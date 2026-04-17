// Per-position coverage and TcellExTRECT-style log-ratio analysis.
//
// The T-cell fraction follows the estimator described in Bentham et al.
// (Nature 2021): within the TCRA locus a VDJ-rearranged cell loses one copy
// of a focal V-J segment, producing a dip in the log2 read-depth ratio
// relative to the surrounding baseline.  The fraction of rearranged cells is
//     f ≈ 1 − 2^logR
// where logR is the median log2(focal / baseline) depth ratio.  The same
// estimator generalises to other V(D)J loci (TRB/TRG/IGH/IGK/IGL).

use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::errors::Error as HtslibError;

pub struct TcellResult {
    pub baseline_median: f64,
    pub focal_median: f64,
    pub log2_ratio: f64,
    pub tcell_fraction: f64,
    pub tcell_fraction_lwr: f64,
    pub tcell_fraction_upr: f64,
    pub focal_n: u64,
    pub baseline_n: u64,
}

pub fn region_coverage(
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    min_mapq: u8,
) -> Result<Vec<(i64, u32)>, HtslibError> {
    let mut reader = IndexedReader::from_path(bam_path)?;
    let tid = reader
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| HtslibError::Fetch)?;
    reader.fetch((tid, start as u64, end as u64))?;

    let mut cov = Vec::new();
    let mut pileup = reader.pileup();
    pileup.set_max_depth(1_000_000);

    for p in pileup {
        let pileup = p?;
        let pos = pileup.pos() as i64;
        if pos < start || pos >= end {
            continue;
        }
        let depth = pileup
            .alignments()
            .filter(|a| {
                let r = a.record();
                !r.is_duplicate()
                    && !r.is_secondary()
                    && !r.is_supplementary()
                    && !r.is_quality_check_failed()
                    && r.mapq() >= min_mapq
                    && !a.is_del()
                    && !a.is_refskip()
            })
            .count() as u32;
        cov.push((pos, depth));
    }
    Ok(cov)
}

pub fn tcell_fraction(
    bam_path: &str,
    chrom: &str,
    focal: (i64, i64),
    left_baseline: (i64, i64),
    right_baseline: (i64, i64),
    min_mapq: u8,
    min_cov: u32,
    exons: &[(i64, i64)],
) -> Result<TcellResult, HtslibError> {
    let lo = left_baseline.0.min(focal.0);
    let hi = right_baseline.1.max(focal.1);
    let cov = region_coverage(bam_path, chrom, lo, hi, min_mapq)?;
    Ok(fraction_from_positions(
        &cov, focal, left_baseline, right_baseline, min_cov, exons,
    ))
}

/// Pure compute on an already-extracted coverage track — the form
/// TcellExTRECT's R `runTcellExTRECT` takes.  Useful for cross-validating
/// against their reference `cov_example` fixture without touching a BAM.
pub fn fraction_from_positions(
    cov: &[(i64, u32)],
    focal: (i64, i64),
    left_baseline: (i64, i64),
    right_baseline: (i64, i64),
    min_cov: u32,
    exons: &[(i64, i64)],
) -> TcellResult {
    let in_range = |pos: i64, (s, e): (i64, i64)| pos >= s && pos < e;
    // If a capture-kit exon list is supplied, restrict coverage to positions
    // that fall inside any exon interval.  This is required for WES BAMs: off-
    // target bases are near-zero and bait-edge positions bias median depth
    // asymmetrically between focal and baseline windows unless both are
    // clipped to the same captured set.  An empty `exons` slice preserves
    // legacy whole-window behaviour (appropriate for WGS).
    let in_exons = |pos: i64| {
        if exons.is_empty() {
            return true;
        }
        exons.iter().any(|&(s, e)| pos >= s && pos < e)
    };

    let mut baseline: Vec<f64> = cov
        .iter()
        .filter(|(p, d)| {
            *d >= min_cov
                && in_exons(*p)
                && (in_range(*p, left_baseline) || in_range(*p, right_baseline))
        })
        .map(|(_, d)| *d as f64)
        .collect();
    let mut focal_vals: Vec<f64> = cov
        .iter()
        .filter(|(p, d)| *d >= min_cov && in_exons(*p) && in_range(*p, focal))
        .map(|(_, d)| *d as f64)
        .collect();

    if baseline.is_empty() || focal_vals.is_empty() {
        return TcellResult {
            baseline_median: f64::NAN,
            focal_median: f64::NAN,
            log2_ratio: f64::NAN,
            tcell_fraction: f64::NAN,
            tcell_fraction_lwr: f64::NAN,
            tcell_fraction_upr: f64::NAN,
            focal_n: focal_vals.len() as u64,
            baseline_n: baseline.len() as u64,
        };
    }

    let baseline_median = median(&mut baseline);
    let focal_median = median(&mut focal_vals);
    let log2_ratio = (focal_median / baseline_median).log2();

    // Per-position log2 ratios, used for a non-parametric 95% CI
    // (2.5% / 97.5% percentiles, bounded to [0, 1] after transform).
    let mut focal_log_ratios: Vec<f64> = focal_vals
        .iter()
        .map(|v| (*v / baseline_median).log2())
        .filter(|x| x.is_finite())
        .collect();
    focal_log_ratios.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let (lwr_lr, upr_lr) = if focal_log_ratios.is_empty() {
        (f64::NAN, f64::NAN)
    } else {
        let n = focal_log_ratios.len();
        let lo_idx = ((n as f64) * 0.025).floor() as usize;
        let hi_idx = (((n as f64) * 0.975).ceil() as usize).saturating_sub(1).min(n - 1);
        (focal_log_ratios[lo_idx], focal_log_ratios[hi_idx])
    };

    let fraction = clamp01(1.0 - 2f64.powf(log2_ratio));
    // Note the inversion: a *lower* logR corresponds to a *higher* fraction.
    let fraction_upr = clamp01(1.0 - 2f64.powf(lwr_lr));
    let fraction_lwr = clamp01(1.0 - 2f64.powf(upr_lr));

    TcellResult {
        baseline_median,
        focal_median,
        log2_ratio,
        tcell_fraction: fraction,
        tcell_fraction_lwr: fraction_lwr,
        tcell_fraction_upr: fraction_upr,
        focal_n: focal_vals.len() as u64,
        baseline_n: baseline.len() as u64,
    }
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = values.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) * 0.5
    }
}

fn clamp01(x: f64) -> f64 {
    if !x.is_finite() {
        x
    } else if x < 0.0 {
        0.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_handles_empty_odd_even() {
        assert!(median(&mut []).is_nan());
        assert_eq!(median(&mut [1.0, 2.0, 3.0]), 2.0);
        assert_eq!(median(&mut [1.0, 2.0, 3.0, 4.0]), 2.5);
    }

    #[test]
    fn fraction_clamping() {
        assert_eq!(clamp01(-0.5), 0.0);
        assert_eq!(clamp01(1.5), 1.0);
        assert_eq!(clamp01(0.3), 0.3);
    }
}
