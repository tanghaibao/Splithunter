// Minimal FFI surface over BWA-MEM.
//
// bwa-sys compiles lh3/bwa's C source into a static libbwa.a linked into this
// crate.  We declare the C symbols we actually need here; we deliberately keep
// opaque structs opaque (via zero-sized `_unused` markers) wherever we don't
// need to inspect their fields — this keeps ABI churn isolated to BWA's C
// definitions of `bwaidx_t`, `mem_alnreg_t`, and `mem_aln_t`.

#![allow(non_camel_case_types)]

use std::os::raw::{c_char, c_int, c_uchar, c_uint, c_ulong};

// --- Opaque handles ---------------------------------------------------------
#[repr(C)]
pub struct bwt_t {
    _unused: [u8; 0],
}
#[repr(C)]
pub struct bntseq_t {
    _unused: [u8; 0],
}
#[repr(C)]
pub struct mem_opt_t {
    // BWA defines many fields.  We never mutate them; we only treat this as
    // an opaque handle produced by mem_opt_init().
    _unused: [u8; 0],
}

/// `bwaidx_t` from bwa/bwa.h.  Fields are all we need.
#[repr(C)]
pub struct bwaidx_t {
    pub bwt: *mut bwt_t,
    pub bns: *mut bntseq_t,
    pub pac: *mut c_uchar,
    pub is_shm: c_int,
    pub l_mem: i64,
    pub mem: *mut c_uchar,
}

/// `mem_alnreg_t` from bwa/bwamem.h.  Field order matches BWA 0.7.17 / master.
#[repr(C)]
#[derive(Copy, Clone)]
pub struct mem_alnreg_t {
    pub rb: i64,
    pub re: i64,
    pub qb: c_int,
    pub qe: c_int,
    pub rid: c_int,
    pub score: c_int,
    pub truesc: c_int,
    pub sub: c_int,
    pub alt_sc: c_int,
    pub csub: c_int,
    pub sub_n: c_int,
    pub w: c_int,
    pub seedcov: c_int,
    pub secondary: c_int,
    pub secondary_all: c_int,
    pub seedlen0: c_int,
    /// bit-packed: n_comp:30, is_alt:2 — accessed via helpers.
    pub n_comp_is_alt: c_int,
    pub frac_rep: f32,
    pub hash: u64,
}

/// Kvec-style vector returned by `mem_align1`.
#[repr(C)]
pub struct mem_alnreg_v {
    pub n: c_ulong,
    pub m: c_ulong,
    pub a: *mut mem_alnreg_t,
}

/// `mem_aln_t` from bwa/bwamem.h — the convenience struct used for reporting.
#[repr(C)]
pub struct mem_aln_t {
    pub pos: i64,
    pub rid: c_int,
    pub flag: c_int,
    /// Bit-packed: is_rev:1, is_alt:1, mapq:8, NM:22.
    pub packed: c_uint,
    pub n_cigar: c_int,
    pub cigar: *mut c_uint,
    pub xa: *mut c_char,
    pub score: c_int,
    pub sub: c_int,
    pub alt_sc: c_int,
}

impl mem_aln_t {
    #[inline]
    pub fn mapq(&self) -> u32 {
        (self.packed >> 2) & 0xff
    }
    #[inline]
    pub fn is_rev(&self) -> bool {
        (self.packed & 1) != 0
    }
}

// --- Extern symbols ---------------------------------------------------------
//
// See bwa/bwa.h and bwa/bwamem.h for the authoritative prototypes.
pub const BWA_IDX_ALL: c_int = 0x7;
pub const BWTALGO_AUTO: c_int = 0;

extern "C" {
    /// Build a BWA-MEM index at `prefix` from the given FASTA.  Returns 0 on
    /// success.  Writes `prefix.{bwt,pac,ann,amb,sa}` to disk.
    pub fn bwa_idx_build(
        fa: *const c_char,
        prefix: *const c_char,
        algo_type: c_int,
        block_size: c_int,
    ) -> c_int;

    /// Load a pre-built BWA index from disk.  Caller owns the returned pointer
    /// and must free with `bwa_idx_destroy`.
    pub fn bwa_idx_load_from_disk(hint: *const c_char, which: c_int) -> *mut bwaidx_t;

    /// Free a `bwaidx_t` returned by `bwa_idx_load_from_disk`.
    pub fn bwa_idx_destroy(idx: *mut bwaidx_t);

    /// Allocate a default `mem_opt_t`.  Caller owns the returned pointer and
    /// must free with `libc::free`.
    pub fn mem_opt_init() -> *mut mem_opt_t;

    /// Align a single short read against the loaded index, returning a kvec of
    /// local alignment regions.  The caller must free `regs.a` with `libc::free`
    /// after consuming the results.
    pub fn mem_align1(
        opt: *const mem_opt_t,
        bwt: *const bwt_t,
        bns: *const bntseq_t,
        pac: *const c_uchar,
        l_seq: c_int,
        seq: *const c_char,
    ) -> mem_alnreg_v;

    /// Convert an alignment region into a reportable `mem_aln_t` (includes
    /// CIGAR, mapq, NM, etc.).  Caller must free `aln.cigar` and `aln.XA` with
    /// `libc::free`.
    pub fn mem_reg2aln(
        opt: *const mem_opt_t,
        bns: *const bntseq_t,
        pac: *const c_uchar,
        l_query: c_int,
        query: *const c_char,
        ar: *const mem_alnreg_t,
    ) -> mem_aln_t;
}
