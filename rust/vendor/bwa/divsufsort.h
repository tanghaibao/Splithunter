/*
 * Stub divsufsort.h — bwtindex.c includes this unconditionally but on small
 * references (< 50 MB, which covers every TCR/IG locus) BWA takes the
 * `is_bwt` / `rope` code path and never calls into libdivsufsort.
 * The empty include lets the file compile without adding a 3rd-party dep.
 */
#ifndef DIVSUFSORT_H
#define DIVSUFSORT_H
#endif
