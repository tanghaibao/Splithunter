// Rust core for Splithunter.
//
// Two public entry points are exposed to Python:
//   analyze_locus  — split-read / split-pair discovery in a BAM region.
//   coverage_logr  — per-exon coverage log-ratio + TcellExTRECT-style fraction.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

mod entropy;
mod coverage;
mod split;
mod bwa_ffi;
mod realign;

/// Analyze a genomic region of a BAM for split-read (SR) and split-pair (SP)
/// evidence.  The returned dict mirrors the JSON schema emitted by the legacy
/// C++ Splithunter binary:
///     <name>.SR-TOTAL, .SR-SIGNAL, .SR-PPM, .SR-DETAILS
///     <name>.SP-TOTAL, .SP-SIGNAL, .SP-PPM, .SP-DETAILS
///
/// `reference_fasta` must be an indexed (FAI) FASTA; the locus slice it
/// yields is used to build a per-call in-memory BWA-MEM index for re-aligning
/// soft-clipped reads, matching the C++ tool's fidelity.
#[pyfunction]
#[pyo3(signature = (bam_path, name, chrom, start, end, reference_fasta, pad=30, indel=10_000, min_entropy=50.0))]
#[allow(clippy::too_many_arguments)]
fn analyze_locus<'py>(
    py: Python<'py>,
    bam_path: &str,
    name: &str,
    chrom: &str,
    start: i64,
    end: i64,
    reference_fasta: &str,
    pad: i64,
    indel: i64,
    min_entropy: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let realigner = realign::Realigner::from_reference(reference_fasta, chrom, start, end, name)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
    let summary = split::analyze(bam_path, chrom, start, end, pad, indel, min_entropy, &realigner)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    let sr_ppm = summary.sr_ppm();
    let sp_ppm = summary.sp_ppm();
    let dict = PyDict::new_bound(py);
    dict.set_item(format!("{name}.SR-TOTAL"), summary.sr_total)?;
    dict.set_item(format!("{name}.SR-SIGNAL"), summary.sr_signal)?;
    dict.set_item(format!("{name}.SR-PPM"), sr_ppm)?;
    dict.set_item(format!("{name}.SR-DETAILS"), summary.sr_details)?;
    dict.set_item(format!("{name}.SP-TOTAL"), summary.sp_total)?;
    dict.set_item(format!("{name}.SP-SIGNAL"), summary.sp_signal)?;
    dict.set_item(format!("{name}.SP-PPM"), sp_ppm)?;
    dict.set_item(format!("{name}.SP-DETAILS"), summary.sp_details)?;
    Ok(dict)
}

/// Per-position coverage across a region; returned as a list of (pos, depth)
/// tuples.  Filters duplicate and secondary alignments by default.
#[pyfunction]
#[pyo3(signature = (bam_path, chrom, start, end, min_mapq=1))]
fn region_coverage<'py>(
    py: Python<'py>,
    bam_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    min_mapq: u8,
) -> PyResult<Bound<'py, PyList>> {
    let cov = coverage::region_coverage(bam_path, chrom, start, end, min_mapq)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
    let list = PyList::empty_bound(py);
    for (pos, depth) in cov {
        list.append((pos, depth))?;
    }
    Ok(list)
}

/// TcellExTRECT-style log-ratio estimate.
///
/// Given a focal (V–J) region flanked by two baseline windows, compute the
/// median log2(focal / baseline) coverage ratio and derive a T-cell fraction
/// estimate (1 − 2^logR, clamped to [0, 1]).
#[pyfunction]
#[pyo3(signature = (bam_path, chrom, focal_start, focal_end, left_start, left_end, right_start, right_end, min_mapq=1, min_cov=1, target_mask=None))]
#[allow(clippy::too_many_arguments)]
fn tcell_fraction<'py>(
    py: Python<'py>,
    bam_path: &str,
    chrom: &str,
    focal_start: i64,
    focal_end: i64,
    left_start: i64,
    left_end: i64,
    right_start: i64,
    right_end: i64,
    min_mapq: u8,
    min_cov: u32,
    target_mask: Option<Vec<(i64, i64)>>,
) -> PyResult<Bound<'py, PyDict>> {
    let target_mask = target_mask.unwrap_or_default();
    let result = coverage::tcell_fraction(
        bam_path,
        chrom,
        (focal_start, focal_end),
        (left_start, left_end),
        (right_start, right_end),
        min_mapq,
        min_cov,
        &target_mask,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    let dict = PyDict::new_bound(py);
    dict.set_item("baseline_median", result.baseline_median)?;
    dict.set_item("focal_median", result.focal_median)?;
    dict.set_item("log2_ratio", result.log2_ratio)?;
    dict.set_item("tcell_fraction", result.tcell_fraction)?;
    dict.set_item("tcell_fraction_lwr", result.tcell_fraction_lwr)?;
    dict.set_item("tcell_fraction_upr", result.tcell_fraction_upr)?;
    dict.set_item("focal_positions", result.focal_n)?;
    dict.set_item("baseline_positions", result.baseline_n)?;
    Ok(dict)
}

/// Pure-compute variant of `tcell_fraction` — takes pre-computed coverage as
/// parallel `positions` + `depths` lists.  Matches the R TcellExTRECT signature
/// where the first argument is a coverage DataFrame.
#[pyfunction]
#[pyo3(signature = (positions, depths, focal_start, focal_end, left_start, left_end, right_start, right_end, min_cov=1, target_mask=None))]
#[allow(clippy::too_many_arguments)]
fn fraction_from_coverage<'py>(
    py: Python<'py>,
    positions: Vec<i64>,
    depths: Vec<u32>,
    focal_start: i64,
    focal_end: i64,
    left_start: i64,
    left_end: i64,
    right_start: i64,
    right_end: i64,
    min_cov: u32,
    target_mask: Option<Vec<(i64, i64)>>,
) -> PyResult<Bound<'py, PyDict>> {
    if positions.len() != depths.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "positions and depths must have the same length",
        ));
    }
    let target_mask = target_mask.unwrap_or_default();
    let cov: Vec<(i64, u32)> = positions.into_iter().zip(depths.into_iter()).collect();
    let result = coverage::fraction_from_positions(
        &cov,
        (focal_start, focal_end),
        (left_start, left_end),
        (right_start, right_end),
        min_cov,
        &target_mask,
    );
    let dict = PyDict::new_bound(py);
    dict.set_item("baseline_median", result.baseline_median)?;
    dict.set_item("focal_median", result.focal_median)?;
    dict.set_item("log2_ratio", result.log2_ratio)?;
    dict.set_item("tcell_fraction", result.tcell_fraction)?;
    dict.set_item("tcell_fraction_lwr", result.tcell_fraction_lwr)?;
    dict.set_item("tcell_fraction_upr", result.tcell_fraction_upr)?;
    dict.set_item("focal_positions", result.focal_n)?;
    dict.set_item("baseline_positions", result.baseline_n)?;
    Ok(dict)
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(analyze_locus, m)?)?;
    m.add_function(wrap_pyfunction!(region_coverage, m)?)?;
    m.add_function(wrap_pyfunction!(tcell_fraction, m)?)?;
    m.add_function(wrap_pyfunction!(fraction_from_coverage, m)?)?;
    Ok(())
}
