// Split-read (SR) and split-pair (SP) detection at a V(D)J locus.
//
// Faithful port of src/Splithunter.cpp: for each non-duplicate read that
// carries a significant soft-clip we re-align the full read against an
// in-memory BWA-MEM index of the locus, split at the query clip junction,
// re-align the two halves independently, then apply the original
// AS / distance / entropy filters.  Split pairs are detected by caching
// paired reads in the target region and checking inter-mate distance.

use std::collections::HashMap;
use std::str;

use rust_htslib::bam::record::{Cigar, Record};
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::errors::Error as HtslibError;

use crate::entropy::trinucleotide_entropy;
use crate::realign::{Hit, Realigner};

#[derive(Default)]
pub struct Summary {
    pub sr_total: u64,
    pub sr_signal: u64,
    pub sr_details: String,
    pub sp_total: u64,
    pub sp_signal: u64,
    pub sp_details: String,
}

impl Summary {
    pub fn sr_ppm(&self) -> f64 {
        if self.sr_total == 0 { 0.0 } else { self.sr_signal as f64 * 1e6 / self.sr_total as f64 }
    }
    pub fn sp_ppm(&self) -> f64 {
        if self.sp_total == 0 { 0.0 } else { self.sp_signal as f64 * 1e6 / self.sp_total as f64 }
    }
}

struct PairedEndpoint {
    chrom: String,
    start: i64,
    end: i64,
    strand: char,
    seq: Vec<u8>,
}

pub fn analyze(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    pad: i64,
    indel: i64,
    min_entropy: f64,
    realigner: &Realigner,
) -> Result<Summary, HtslibError> {
    let mut reader = IndexedReader::from_path(bam_path)?;
    let tid = reader
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| HtslibError::Fetch)?;
    let header_names: Vec<String> = (0..reader.header().target_count())
        .map(|i| String::from_utf8_lossy(reader.header().tid2name(i)).into_owned())
        .collect();
    reader.fetch((tid, start as u64, end as u64))?;

    let mut summary = Summary::default();
    let mut pair_cache: HashMap<String, Vec<PairedEndpoint>> = HashMap::new();
    let mut rec = Record::new();

    while let Some(result) = reader.read(&mut rec) {
        result?;
        if rec.is_duplicate() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }

        summary.sr_total += 1;

        let qname = str::from_utf8(rec.qname()).unwrap_or("").to_string();
        let read_len = rec.seq_len() as i64;
        let clip = soft_clip_len(&rec);

        // ---- SP cache: paired read that is "mostly aligned" ---------------
        //
        // Matches the C++ condition readScore >= readLen - PAD, where score
        // is the number of aligned bases.  Aligned == readLen - soft_clip.
        let aligned_primary = read_len - clip;
        if rec.is_paired()
            && !rec.is_unmapped()
            && !rec.is_mate_unmapped()
            && aligned_primary >= read_len - pad
        {
            let start_ref = rec.pos();
            let end_ref = rec.cigar().end_pos();
            let strand = if rec.is_reverse() { '-' } else { '+' };
            let seq = rec.seq().as_bytes();
            let chrom_name = header_names
                .get(rec.tid() as usize)
                .cloned()
                .unwrap_or_else(|| chrom.to_string());
            pair_cache.entry(qname.clone()).or_default().push(PairedEndpoint {
                chrom: chrom_name,
                start: start_ref,
                end: end_ref,
                strand,
                seq,
            });
        }

        // ---- SR detection -------------------------------------------------
        //
        // Only reads with a significant clip are re-aligned.  The C++ tool
        // checked `NumClip() < PAD || > readLen - PAD` — reject both
        // untrimmed reads (clip too small) and reads that are almost
        // entirely clipped (nothing to align against the locus).
        if clip < pad || clip > read_len - pad {
            continue;
        }

        // Full-read realignment against the locus.
        let read_seq = rec.seq().as_bytes();
        let full_hits = realigner.align(&read_seq);
        let Some(primary) = full_hits.first() else {
            continue;
        };

        // The realigned primary must still carry a significant clip on one
        // side — otherwise the original clip was just noise and the read is
        // not chimeric against this locus.
        let realigned_clip_left = primary.query_start;
        let realigned_clip_right = read_len - primary.query_end;
        let (split_idx, left_part, right_part) = if realigned_clip_left > pad
            && realigned_clip_left < read_len - pad
        {
            // Clipped at the 5' end: left = the clipped prefix, right = the
            // rest starting at the realigned query start.
            let split = primary.query_start as usize;
            (split, &read_seq[..split], &read_seq[split..])
        } else if realigned_clip_right > pad && realigned_clip_right < read_len - pad {
            // Clipped at the 3' end.
            let split = primary.query_end as usize;
            (split, &read_seq[..split], &read_seq[split..])
        } else {
            continue;
        };
        let _ = split_idx; // retained for future diagnostics

        // Re-align each half independently.
        let left_hits = realigner.align(left_part);
        let right_hits = realigner.align(right_part);
        let (Some(l), Some(r)) = (left_hits.first(), right_hits.first()) else {
            continue;
        };

        // SR Condition 1: each part is significant (AS score proxy).
        if l.as_score < pad || r.as_score < pad {
            continue;
        }
        // SR Condition 2: combined score exceeds readLen - PAD/2.
        if l.as_score + r.as_score < read_len - pad / 2 {
            continue;
        }
        // SR Condition 3: distinct regions (V(D)J rearrangement signal).
        let dist = (r.ref_start - l.ref_start).abs();
        if dist < indel {
            continue;
        }
        // SR Condition 4: sequence complexity — each half separately.
        if trinucleotide_entropy(left_part) < min_entropy
            || trinucleotide_entropy(right_part) < min_entropy
        {
            continue;
        }

        summary.sr_signal += 1;
        write_hit_detail(&mut summary.sr_details, l, r);
    }

    // ---- SP pair scan ---------------------------------------------------
    for (_, mut mates) in pair_cache.into_iter() {
        if mates.len() != 2 {
            continue;
        }
        summary.sp_total += 1;
        mates.sort_by_key(|m| m.start);
        let (left, right) = (&mates[0], &mates[1]);
        let dist = (right.start - left.start).abs();
        if dist < indel {
            continue;
        }
        if trinucleotide_entropy(&left.seq) < min_entropy
            || trinucleotide_entropy(&right.seq) < min_entropy
        {
            continue;
        }
        summary.sp_signal += 1;
        append_detail(
            &mut summary.sp_details,
            &left.chrom, left.start, left.end, left.strand,
            &right.chrom, right.start, right.end, right.strand,
        );
    }

    Ok(summary)
}

fn soft_clip_len(rec: &Record) -> i64 {
    let mut clip = 0_i64;
    for op in rec.cigar().iter() {
        if let Cigar::SoftClip(n) = op {
            clip += *n as i64;
        }
    }
    clip
}

fn write_hit_detail(details: &mut String, l: &Hit, r: &Hit) {
    append_detail(
        details,
        &l.contig_name, l.ref_start, l.ref_end, l.strand,
        &r.contig_name, r.ref_start, r.ref_end, r.strand,
    );
}

#[allow(clippy::too_many_arguments)]
fn append_detail(
    details: &mut String,
    left_chrom: &str, left_start: i64, left_end: i64, left_strand: char,
    right_chrom: &str, right_start: i64, right_end: i64, right_strand: char,
) {
    use std::fmt::Write as _;
    let _ = write!(
        details,
        "{}:{}-{}({})|{}:{}-{}({});",
        left_chrom, left_start, left_end, left_strand,
        right_chrom, right_start, right_end, right_strand,
    );
}
