// Pure-FFI wrapper around BWA-MEM: extract a locus slice from a user-supplied
// reference FASTA, build a per-locus BWA index in a tempdir, load it, and
// expose a single-read alignment helper.
//
// The port of Splithunter.cpp relies on three BWA behaviours that minimap2 does
// not faithfully reproduce for short reads: (1) the BWA seed-and-extend with
// its specific soft-clip behaviour, (2) the AS scoring scale used in the
// original `AS >= readLen - PAD/2` filters, and (3) the chimeric supplementary
// split that the original tool drove by re-aligning the soft-clipped tail with
// the same BWA index.  All three come from libbwa directly.

use std::ffi::CString;
use std::os::raw::c_char;
use std::path::{Path, PathBuf};
use std::ptr;

use rust_htslib::faidx::Reader as FaiReader;
use tempfile::TempDir;

use crate::bwa_ffi::*;

#[derive(Debug)]
pub enum RealignError {
    Faidx(String),
    Io(String),
    BuildIndex(i32),
    LoadIndex(String),
}

impl std::fmt::Display for RealignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RealignError::Faidx(s) => write!(f, "faidx: {s}"),
            RealignError::Io(s) => write!(f, "io: {s}"),
            RealignError::BuildIndex(c) => write!(f, "bwa_idx_build returned {c}"),
            RealignError::LoadIndex(s) => write!(f, "bwa_idx_load_from_disk: {s}"),
        }
    }
}
impl std::error::Error for RealignError {}

/// A loaded per-locus BWA-MEM index.
///
/// The backing tempdir is kept alive for the lifetime of the aligner so the
/// `.bwt`/`.pac`/... files remain on disk until we're done with them.
pub struct Realigner {
    idx: *mut bwaidx_t,
    opt: *mut mem_opt_t,
    contig_name: String,
    contig_offset: i64,
    _tmp: TempDir,
}

// BWA's index and option structs are read-only after construction and safe to
// share across threads — we only call const-pointer FFI entry points on them.
unsafe impl Send for Realigner {}
unsafe impl Sync for Realigner {}

impl Drop for Realigner {
    fn drop(&mut self) {
        unsafe {
            if !self.idx.is_null() {
                bwa_idx_destroy(self.idx);
            }
            if !self.opt.is_null() {
                libc::free(self.opt as *mut libc::c_void);
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct Hit {
    pub contig_name: String,
    pub ref_start: i64,
    pub ref_end: i64,
    pub query_start: i64,
    pub query_end: i64,
    pub aligned_bases: i64,
    pub as_score: i64,
    pub strand: char,
    pub mapq: u32,
}

impl Realigner {
    /// Extract `chrom:start-end` from `reference_fasta`, write a single-record
    /// FASTA to a tempdir, run `bwa_idx_build`, then `bwa_idx_load_from_disk`.
    pub fn from_reference(
        reference_fasta: &str,
        chrom: &str,
        start: i64,
        end: i64,
        locus_name: &str,
    ) -> Result<Self, RealignError> {
        let reader = FaiReader::from_path(reference_fasta)
            .map_err(|e| RealignError::Faidx(format!("open {reference_fasta}: {e}")))?;
        // fetch_seq_string takes inclusive end.
        let seq = reader
            .fetch_seq_string(chrom, start as usize, (end - 1) as usize)
            .map_err(|e| RealignError::Faidx(format!("fetch {chrom}:{start}-{end}: {e}")))?;
        if seq.is_empty() {
            return Err(RealignError::Faidx(format!(
                "empty slice for {chrom}:{start}-{end}"
            )));
        }

        let tmp = tempfile::Builder::new()
            .prefix("splithunter-bwa-")
            .tempdir()
            .map_err(|e| RealignError::Io(e.to_string()))?;
        let fasta_path = tmp.path().join(format!("{locus_name}.fa"));
        let record = format!(">{locus_name}\n{seq}\n");
        std::fs::write(&fasta_path, record).map_err(|e| RealignError::Io(e.to_string()))?;

        let prefix_path: PathBuf = fasta_path.clone();
        build_index(&fasta_path, &prefix_path)?;
        let idx = load_index(&prefix_path)?;
        let opt = unsafe { mem_opt_init() };
        if opt.is_null() {
            unsafe { bwa_idx_destroy(idx) };
            return Err(RealignError::LoadIndex("mem_opt_init() returned null".into()));
        }

        Ok(Self {
            idx,
            opt,
            contig_name: locus_name.to_string(),
            contig_offset: start,
            _tmp: tmp,
        })
    }

    /// Align a short read against the locus index.  Returns hits sorted by
    /// alignment score descending (BWA's default).
    pub fn align(&self, read: &[u8]) -> Vec<Hit> {
        if read.is_empty() || self.idx.is_null() {
            return Vec::new();
        }
        let idx = unsafe { &*self.idx };

        let c_seq = match CString::new(read) {
            Ok(v) => v,
            Err(_) => return Vec::new(),
        };

        let regs = unsafe {
            mem_align1(
                self.opt,
                idx.bwt,
                idx.bns,
                idx.pac,
                read.len() as i32,
                c_seq.as_ptr(),
            )
        };
        if regs.a.is_null() || regs.n == 0 {
            if !regs.a.is_null() {
                unsafe { libc::free(regs.a as *mut libc::c_void) };
            }
            return Vec::new();
        }

        let mut hits = Vec::with_capacity(regs.n as usize);
        for i in 0..regs.n {
            let reg = unsafe { *regs.a.add(i as usize) };
            // Skip secondary alignments and low-quality chains.
            if reg.secondary >= 0 {
                continue;
            }
            let aln = unsafe {
                mem_reg2aln(
                    self.opt,
                    idx.bns,
                    idx.pac,
                    read.len() as i32,
                    c_seq.as_ptr(),
                    &reg as *const _,
                )
            };
            if aln.rid < 0 {
                free_aln(&aln);
                continue;
            }
            // `reg.rb/re` point into BWA's packed forward+reverse reference
            // and are NOT real coordinates for reverse-strand hits; only
            // `aln.pos` is canonical.  Use aln.pos as the forward-strand
            // reference start and preserve the reference span (re - rb).
            let ref_span = (reg.re - reg.rb) as i64;
            hits.push(Hit {
                contig_name: self.contig_name.clone(),
                ref_start: self.contig_offset + aln.pos,
                ref_end: self.contig_offset + aln.pos + ref_span,
                query_start: reg.qb as i64,
                query_end: reg.qe as i64,
                aligned_bases: (reg.qe - reg.qb) as i64,
                as_score: reg.score as i64,
                strand: if aln.is_rev() { '-' } else { '+' },
                mapq: aln.mapq(),
            });
            free_aln(&aln);
        }
        unsafe { libc::free(regs.a as *mut libc::c_void) };
        hits
    }
}

fn build_index(fa: &Path, prefix: &Path) -> Result<(), RealignError> {
    let fa_c = path_to_cstring(fa)?;
    let prefix_c = path_to_cstring(prefix)?;
    let rc = unsafe {
        bwa_idx_build(
            fa_c.as_ptr(),
            prefix_c.as_ptr(),
            BWTALGO_AUTO,
            10_000_000,
        )
    };
    if rc != 0 {
        return Err(RealignError::BuildIndex(rc));
    }
    Ok(())
}

fn load_index(prefix: &Path) -> Result<*mut bwaidx_t, RealignError> {
    let prefix_c = path_to_cstring(prefix)?;
    let idx = unsafe { bwa_idx_load_from_disk(prefix_c.as_ptr(), BWA_IDX_ALL) };
    if idx.is_null() {
        return Err(RealignError::LoadIndex(format!(
            "bwa_idx_load_from_disk({}) returned null",
            prefix.display()
        )));
    }
    Ok(idx)
}

fn path_to_cstring(p: &Path) -> Result<CString, RealignError> {
    let s = p
        .to_str()
        .ok_or_else(|| RealignError::Io(format!("non-utf8 path {}", p.display())))?;
    CString::new(s).map_err(|e| RealignError::Io(e.to_string()))
}

fn free_aln(aln: &mem_aln_t) {
    unsafe {
        if !aln.cigar.is_null() {
            libc::free(aln.cigar as *mut libc::c_void);
        }
        if !aln.xa.is_null() {
            libc::free(aln.xa as *mut libc::c_void);
        }
    }
    // Silence "unused field" warnings from the compiler for fields we keep
    // around for ABI completeness only.
    let _ = aln as *const mem_aln_t as *const c_char;
    let _ = ptr::null::<u8>();
}

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end sanity check: build a BWA-MEM index from a synthetic 5 kb
    /// reference, align a read drawn from position 1000, and verify the
    /// realigner returns a hit whose reference coordinate lands in the
    /// correct window.  This proves the entire vendored-BWA + FFI chain is
    /// wired up correctly on the host architecture (SSE2/NEON included).
    #[test]
    fn realigner_roundtrip_on_synthetic_reference() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("toy.fa");
        // Deterministic high-entropy sequence via xorshift — short cycles
        // in a toy reference would give BWA equally-scoring hits and mask
        // position correctness.
        let mut seq = String::with_capacity(5000);
        let alphabet = b"ACGT";
        let mut state: u64 = 0x9E3779B97F4A7C15;
        for _ in 0..5000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            seq.push(alphabet[(state as usize) & 3] as char);
        }
        std::fs::write(&fa, format!(">chr_toy\n{seq}\n")).unwrap();
        // Build a FAI so the faidx reader can slice.
        let (out, _) = std::process::Command::new("/bin/sh")
            .arg("-c")
            .arg(format!(
                "which samtools >/dev/null && samtools faidx {} && echo fai",
                fa.display()
            ))
            .output()
            .map(|o| (o.stdout, o.stderr))
            .unwrap_or_default();
        if !String::from_utf8_lossy(&out).contains("fai") {
            // Synthesize a .fai manually: <name>\t<len>\t<offset>\t<linebases>\t<linewidth>
            let header_len = "chr_toy".len() as i64 + 2; // ">chr_toy\n"
            let content = format!("chr_toy\t{}\t{}\t{}\t{}\n", seq.len(), header_len, seq.len(), seq.len() + 1);
            std::fs::write(fa.with_extension("fa.fai"), content).unwrap();
        }

        let r = Realigner::from_reference(
            fa.to_str().unwrap(),
            "chr_toy",
            0,
            seq.len() as i64,
            "chr_toy",
        )
        .expect("index build");

        // Grab a 100 bp read starting at position 1000.
        let read: Vec<u8> = seq.as_bytes()[1000..1100].to_vec();
        let hits = r.align(&read);
        assert!(!hits.is_empty(), "no hit on self-aligned read");
        let h = &hits[0];
        assert!(
            (h.ref_start - 1000).abs() < 20,
            "expected ref_start ≈ 1000, got {}",
            h.ref_start
        );
        assert!(h.aligned_bases >= 95, "query span {}", h.aligned_bases);
        assert_eq!(h.strand, '+');
    }
}
