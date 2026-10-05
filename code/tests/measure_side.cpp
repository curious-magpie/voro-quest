// rvd-measure-side: the numbers the documents quote about the side predicates
// and the filter, measured, so that they can be measured again.
//
// Not a test: it checks nothing and is not run by ctest. It prints
//
//   rates     on the input families of side_support.h, with the suites'
//             seeds: how often the filter decides a random input, how often
//             a lattice input is an exact tie, how often a tie nudged by one
//             ulp reaches the exact path, and how often constructing the
//             vertex in doubles gets that nudged tie wrong;
//   timings   nanoseconds per call of each formula evaluated as double, as
//             Bounded and as DynamicExpansion, and of the whole predicate
//             (filter, then exact when needed), on 20000 random inputs, best
//             of 5 runs. Microbenchmarks: the ratios are what to read, and
//             they depend on the machine and the compiler.
//
// Build it in Release (the default) and run build/rvd-measure-side. Record
// the date, the CPU and the compiler with the numbers.
#include <algorithm>
#include <chrono>
#include <cstdio>

#include "side_support.h"

#include "kernel/orient3d.h"
#include "kernel/side.h"

namespace
{

using Clock = std::chrono::steady_clock;

const std::size_t kSites[4] = {2, 3, 4, 5}, kCorners[4] = {1, 2, 3, 0};

int predicate(int p, const Case &c)
{
    const std::vector<Site> &s = c.s;
    switch (p)
    {
    case 0:
        return side1(s[0].p, s[0].w, s[1].p, s[1].w, c.q[0]);
    case 1:
        return side2(s[0].p, s[0].w, s[1].p, s[1].w, s[2].p, s[2].w, c.q[0],
                     c.q[1]);
    case 2:
        return side3(s[0].p, s[0].w, s[1].p, s[1].w, s[2].p, s[2].w, s[3].p,
                     s[3].w, c.q[0], c.q[1], c.q[2]);
    default:
        return side4(s[0].p, s[0].w, s[1].p, s[1].w, s[2].p, s[2].w, s[3].p,
                     s[3].w, s[4].p, s[4].w);
    }
}

PredicateCounts &counts(int p)
{
    switch (p)
    {
    case 0:
        return side1_counts();
    case 1:
        return side2_counts();
    case 2:
        return side3_counts();
    default:
        return side4_counts();
    }
}

// The predicate's numerator and denominator in the number type T.
template <class T> std::pair<T, T> terms(int p, const Case &c)
{
    const std::vector<Site> &s = c.s;
    switch (p)
    {
    case 0:
        return {bisector_value<T>(s[0].p, s[0].w, s[1].p, s[1].w, c.q[0]),
                T(1.0)};
    case 1:
        return side2_terms<T>(s[0].p, s[0].w, s[1].p, s[1].w, s[2].p, s[2].w,
                              c.q[0], c.q[1]);
    case 2:
        return side3_terms<T>(s[0].p, s[0].w, s[1].p, s[1].w, s[2].p, s[2].w,
                              s[3].p, s[3].w, c.q[0], c.q[1], c.q[2]);
    default:
        return side4_terms<T>(s[0].p, s[0].w, s[1].p, s[1].w, s[2].p, s[2].w,
                              s[3].p, s[3].w, s[4].p, s[4].w);
    }
}

// Something to add up from each result, so the compiler keeps the work.
double use(double x)
{
    return x;
}
double use(const Bounded &b)
{
    return b.value + b.error;
}
double use(const DynamicExpansion &e)
{
    return approximate(e);
}

volatile double sink = 0;

template <class F> double best_ns(F f, std::size_t calls)
{
    double best = 1e300;
    for (int r = 0; r < 5; ++r)
    {
        const auto t0 = Clock::now();
        f();
        const std::chrono::duration<double, std::nano> d = Clock::now() - t0;
        best = std::min(best, d.count() / double(calls));
    }
    return best;
}

// The same draws as check_predicate in side_support.h, with its seed.
void rates(int p)
{
    std::mt19937_64 rng(100 * (p + 1));
    counts(p) = PredicateCounts{};
    long n = 0;
    for (int it = 0; it < 50000; ++it)
    {
        const Case c = make_case(rng, kSites[p], kCorners[p], false);
        if (!oracle(c))
            continue;
        ++n;
        predicate(p, c);
    }
    const double filtered = 100.0 * double(counts(p).filtered) / double(n);

    long lattice = 0, ties = 0;
    for (int it = 0; it < 50000; ++it)
    {
        const Case c = make_case(rng, kSites[p], kCorners[p], true);
        const std::optional<int> o = oracle(c);
        if (!o)
            continue;
        ++lattice;
        ties += *o == 0;
    }

    counts(p) = PredicateCounts{};
    long nudges = 0, naive_wrong = 0;
    for (int it = 0; it < 200000 && nudges < 20000; ++it)
    {
        const Case tie = make_case(rng, kSites[p], kCorners[p], true);
        const std::optional<int> at_tie = oracle(tie);
        if (!at_tie || *at_tie != 0)
            continue;
        const Case c = nudged(tie, rng);
        const std::optional<int> want = oracle(c);
        if (!want || *want == 0)
            continue;
        ++nudges;
        predicate(p, c);
        const std::optional<int> guess = naive(c);
        naive_wrong += !guess || *guess != *want;
    }
    std::printf("side%d  random %ld: filter %.2f%%   lattice %ld: ties %.2f%%"
                "   nudged %ld: exact %.1f%%, doubles wrong %.1f%%\n",
                p + 1, n, filtered, lattice, 100.0 * double(ties) / lattice,
                nudges, 100.0 * double(counts(p).exact) / double(nudges),
                100.0 * double(naive_wrong) / double(nudges));
}

void timings(int p)
{
    std::mt19937_64 rng(900 + p);
    std::vector<Case> cs;
    while (cs.size() < 20000)
    {
        Case c = make_case(rng, kSites[p], kCorners[p], false);
        if (oracle(c))
            cs.push_back(c);
    }
    const std::size_t n = cs.size();
    auto formula = [&](auto zero)
    {
        using T = decltype(zero);
        double s = 0;
        for (const Case &c : cs)
        {
            const std::pair<T, T> t = terms<T>(p, c);
            s += use(t.first) + use(t.second);
        }
        sink = s;
    };
    const double d = best_ns([&] { formula(0.0); }, n);
    const double b = best_ns([&] { formula(Bounded()); }, n);
    const double e = best_ns([&] { formula(DynamicExpansion()); }, n);
    const double whole = best_ns(
        [&]
        {
            long s = 0;
            for (const Case &c : cs)
                s += predicate(p, c);
            sink = double(s);
        },
        n);
    std::printf("side%d  double %6.1f   Bounded %6.1f (x%4.1f)   "
                "DynamicExpansion %7.0f (x%5.0f)   predicate %6.1f   ns\n",
                p + 1, d, b, b / d, e, e / d, whole);
}

void orient3d_timings()
{
    std::mt19937_64 rng(950);
    const std::size_t n = 20000;
    std::vector<double> v(12 * n);
    for (double &x : v)
        x = random_double(rng, -20, 20);
    auto formula = [&](auto zero)
    {
        using T = decltype(zero);
        double s = 0;
        for (std::size_t k = 0; k < n; ++k)
        {
            const double *a = &v[12 * k];
            s += use(orient3d_det(T(a[0]), T(a[1]), T(a[2]), T(a[3]), T(a[4]),
                                  T(a[5]), T(a[6]), T(a[7]), T(a[8]), T(a[9]),
                                  T(a[10]), T(a[11])));
        }
        sink = s;
    };
    const double d = best_ns([&] { formula(0.0); }, n);
    const double b = best_ns([&] { formula(Bounded()); }, n);
    const double whole = best_ns(
        [&]
        {
            long s = 0;
            for (std::size_t k = 0; k < n; ++k)
            {
                const double *a = &v[12 * k];
                s += orient3d(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7],
                              a[8], a[9], a[10], a[11]);
            }
            sink = double(s);
        },
        n);
    std::printf("orient3d double %5.1f   Bounded %6.1f (x%4.1f)   "
                "predicate %6.1f   ns\n",
                d, b, b / d, whole);
}

} // namespace

int main()
{
    std::printf("rates, on the side suites' input families and seeds\n");
    for (int p = 0; p < 4; ++p)
        rates(p);
    std::printf("\ntimings, ns per call, 20000 random inputs, best of 5\n");
    for (int p = 0; p < 4; ++p)
        timings(p);
    orient3d_timings();
    return 0;
}
