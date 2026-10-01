#pragma once
#include <cstdint>

// How each call of a predicate was decided: by the filter, or by the exact
// path (design, section 6.6). The Robustness window will show them, and the
// tests use them now, to check that the filter does its job -- deciding the
// easy cases -- and only that.
//
// One set per predicate, per thread: orient2d_counts() returns this thread's
// own, so counting needs no lock and threads never race on it. Totals are a
// sum over threads, which is exact, being integers. Reset with
// `orient2d_counts() = PredicateCounts{};`.
struct PredicateCounts
{
    uint64_t filtered = 0; // decided by the filter
    uint64_t exact = 0;    // needed the exact path
};
