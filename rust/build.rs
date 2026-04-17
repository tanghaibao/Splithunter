// Compile a minimal subset of BWA-MEM's C source in-tree.
//
// We vendor `bwa/` from lh3/bwa (BSD-3) with two local tweaks:
//   * ksw.c: `#include <emmintrin.h>` rewritten to `"splithunter_simd.h"` so
//     the file compiles under Aarch64 (via sse2neon.h).
//   * divsufsort.h: stubbed; bwtindex.c only takes that code path for
//     references > 50 MB, which we never hit for V(D)J loci.

const SOURCES: &[&str] = &[
    "vendor/bwa/utils.c",
    "vendor/bwa/kthread.c",
    "vendor/bwa/kstring.c",
    "vendor/bwa/ksw.c",
    "vendor/bwa/bwt.c",
    "vendor/bwa/bntseq.c",
    "vendor/bwa/bwa.c",
    "vendor/bwa/bwamem.c",
    "vendor/bwa/bwamem_pair.c",
    "vendor/bwa/bwamem_extra.c",
    "vendor/bwa/malloc_wrap.c",
    "vendor/bwa/bwtindex.c",
    "vendor/bwa/is.c",
    "vendor/bwa/rle.c",
    "vendor/bwa/rope.c",
];

fn main() {
    for f in SOURCES {
        println!("cargo:rerun-if-changed={f}");
    }
    println!("cargo:rerun-if-changed=vendor/bwa/splithunter_simd.h");
    println!("cargo:rerun-if-changed=vendor/bwa/divsufsort.h");

    let mut build = cc::Build::new();
    build
        .files(SOURCES)
        .include("vendor/bwa")
        .warnings(false)
        .extra_warnings(false)
        .flag_if_supported("-fPIC")
        .flag_if_supported("-fcommon")
        .define("HAVE_PTHREAD", None);

    // BWA's libz dep.
    println!("cargo:rustc-link-lib=z");

    build.compile("bwa");
}
