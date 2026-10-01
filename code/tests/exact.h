// The oracle's arithmetic: exact rationals, for checking the kernel.
//
// Every double is a rational number, and GMP's mpq_class holds any rational
// exactly and adds and multiplies without rounding. So "is this expansion
// exactly a + b?" becomes a plain comparison of two mpq_class values -- slow,
// and correct by construction, which is all a reference has to be.
//
// Also here: the random doubles the kernel suites draw, spread over exponents
// on purpose. Uniform doubles in [0, 1) would exercise almost none of what
// makes floating-point addition hard: cancellation, absorption of a tiny term,
// operands whose exponents are far apart.
#pragma once

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>

#include <gmpxx.h>

// The exact value of a finite double. mpq_class's double constructor is exact.
inline mpq_class exact(double x)
{
    return mpq_class(x);
}

// How many bits of `x`'s significand lie between its highest and lowest set
// bit, inclusive: 1 for a power of two, 53 at most, 0 for zero. What "fits in
// 26 bits" means for a Dekker split.
inline int significant_bits(double x)
{
    if (x == 0.0)
        return 0;
    int e = 0;
    const double m = std::frexp(std::fabs(x), &e); // in [0.5, 1)
    uint64_t bits = uint64_t(std::ldexp(m, 53));   // exact: 53 bits
    int trailing = 0;
    while ((bits & 1u) == 0)
    {
        bits >>= 1;
        ++trailing;
    }
    return 53 - trailing;
}

// The exponent of `x`'s highest and lowest set bits: x = 1.5 has its highest
// at 0 (the 1) and its lowest at -1 (the .5). Two doubles are nonoverlapping
// when the lowest set bit of the larger is above the highest of the smaller.
inline int highest_bit(double x)
{
    return std::ilogb(x);
}
inline int lowest_bit(double x)
{
    return std::ilogb(x) - significant_bits(x) + 1;
}

// A double with a random sign, a random full 53-bit significand, and an
// exponent drawn uniformly from [emin, emax].
inline double random_double(std::mt19937_64 &rng, int emin, int emax)
{
    std::uniform_int_distribution<int> exponent(emin, emax);
    const uint64_t significand = (rng() >> 11) | (uint64_t(1) << 52);
    const double x = std::ldexp(double(significand), exponent(rng) - 52);
    return (rng() & 1u) ? -x : x;
}

// `x` printed so that it reads back as the same double, for failure messages.
inline std::string show(double x)
{
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.17g", x);
    return buf;
}
