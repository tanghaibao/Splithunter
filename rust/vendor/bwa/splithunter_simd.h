/*
 * Thin shim BWA's ksw.c includes (we rewrote the `<emmintrin.h>` there to
 * point here).  On x86 we forward to the real SSE2 header; on ARM we hand off
 * to sse2neon, which maps Intel SSE intrinsics onto Arm/Aarch64 NEON.
 */
#ifndef SPLITHUNTER_SIMD_H
#define SPLITHUNTER_SIMD_H

#if defined(__aarch64__) || defined(__arm__) || defined(_M_ARM64)
#  ifndef SSE2NEON_ALLOC_DEFINED
#    define SSE2NEON_ALLOC_DEFINED 1
#  endif
#  include "sse2neon.h"
#else
#  include <emmintrin.h>
#endif

#endif
