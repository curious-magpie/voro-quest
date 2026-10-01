// What the expansion suites share: the oracle's view of an expansion (its
// exact value, and the checks of its form), the failure tally, and random
// expansions built with grow().
//
// Header-only and inline, like test_support.h, so each suite stays one source
// file. Everything here is checked by test_expansion.cpp before any other
// suite relies on it: build() trusts grow(), and grow() is tested first.
#pragma once

#include <cmath>
#include <random>
#include <string>

#include "test_support.h"

#include "exact.h"
#include "kernel/expansion.h"

// --- the oracle's view of an expansion ---

// These take any expansion type -- Expansion<N> or DynamicExpansion --
// through the two things both have: size() and operator[].
template <class E> inline mpq_class value(const E &e)
{
    mpq_class v = 0;
    for (std::size_t i = 0; i < e.size(); ++i)
        v += exact(e[i]);
    return v;
}

// Nonzero, and each component entirely below the next in its binary digits.
// Checking neighbours is enough: nonoverlapping is transitive along a list in
// increasing magnitude.
//
// (Not checked here: that size() is within a fixed capacity. append() asserts
// that on every component, and the tests keep asserts on.)
template <class E> inline bool well_formed(const E &e)
{
    for (std::size_t i = 0; i < e.size(); ++i)
        if (e[i] == 0.0 || !std::isfinite(e[i]))
            return false;
    for (std::size_t i = 1; i < e.size(); ++i)
        if (!(highest_bit(e[i - 1]) < lowest_bit(e[i])))
            return false;
    return true;
}

// Two neighbouring components x < y are *adjacent* (Shewchuk, section 2.3)
// when they do not overlap but touch: x's highest set bit is the one just below
// y's lowest, so that 2x would overlap y.
inline bool touching(double x, double y)
{
    return highest_bit(x) == lowest_bit(y) - 1;
}

inline bool power_of_two(double x)
{
    return significant_bits(x) == 1;
}

// Nonadjacent: well formed, and no two neighbours touch. What grow(), sum()
// and scale() produce under round-to-even (Shewchuk's Theorems 10, 12, 19).
template <class E> inline bool nonadjacent(const E &e)
{
    for (std::size_t i = 1; i < e.size(); ++i)
        if (touching(e[i - 1], e[i]))
            return false;
    return well_formed(e);
}

// Strongly nonoverlapping (Shewchuk, section 2.3): well formed, and where two
// neighbours do touch, both are powers of two, and no component touches both
// of its neighbours. Weaker than nonadjacent. It is what Fast-Expansion-Sum
// needs of its inputs, and what it promises of its output (Theorem 13), so a
// fast_sum result can be fed back into fast_sum.
template <class E> inline bool strongly_nonoverlapping(const E &e)
{
    if (!well_formed(e))
        return false;
    for (std::size_t i = 1; i < e.size(); ++i)
    {
        if (!touching(e[i - 1], e[i]))
            continue;
        if (!power_of_two(e[i - 1]) || !power_of_two(e[i]))
            return false;
        if (i + 1 < e.size() && touching(e[i], e[i + 1]))
            return false;
    }
    return true;
}

template <class E> inline std::string describe(const E &e)
{
    std::string s = "[";
    for (std::size_t i = 0; i < e.size(); ++i)
        s += (i ? ", " : "") + show(e[i]);
    return s + "]";
}

// Failures of one property over many inputs, with the first one kept.
struct Tally
{
    int failures = 0;
    std::string first;

    void add(bool ok, const std::string &example)
    {
        if (!ok && failures++ == 0)
            first = example;
    }
    void report(const char *what) const
    {
        check(failures == 0, what, failures ? first.c_str() : nullptr);
    }
};

// --- random inputs ---

// Exponents kept within +-200, so that scaling by factors within +-100 never
// makes an error term underflow (the range kernel/eft.h promises exactness on).
inline double input(std::mt19937_64 &rng)
{
    return random_double(rng, -200, 200);
}

// An expansion grown from K random doubles, and the exact sum of those doubles.
template <std::size_t K>
Expansion<K> build(std::mt19937_64 &rng, mpq_class &sum)
{
    if constexpr (K == 1)
    {
        const double a = input(rng);
        sum = exact(a);
        return Expansion<1>(a);
    }
    else
    {
        Expansion<K - 1> e = build<K - 1>(rng, sum);
        // Every other term nearly cancels what is there, so that results with
        // fewer components than inputs -- and zero -- turn up.
        double b = input(rng);
        if (K % 2 == 0 && e.size() > 0)
            b = -approximate(e) * (1.0 + std::ldexp(1.0, -40));
        sum += exact(b);
        return grow(e, b);
    }
}
