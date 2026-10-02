// What the side-predicate suites share: the oracle, the inputs, the driver.
//
// side1..side4 (kernel/side.h) each answer one question: at a vertex of a
// clipped piece, is seed i's power distance smaller than seed m's? The
// vertex is given by what defines it -- a domain corner, an edge crossing a
// bisector, a triangle crossing two, or four seeds -- never by coordinates.
// The predicates get the sign from a ratio of determinants of f values
// (design, section 6.2). The oracle here does the obvious thing instead, as
// the design asks (section 6.5): it constructs the vertex exactly, in
// rationals, by solving its defining linear system, and compares the two
// power distances straight from their definition, |x - p|^2 - w. It shares no
// formula with the predicates, which is the point: two derivations from
// different starting points that agree on millions of inputs.
//
// The same construction run in doubles is what a floating-point Voronoi code
// does, and the suites use it to show that their hard inputs are hard: it
// gets some of them wrong.
//
// One description covers all four predicates. A Case holds the seeds, in the
// order i, (j, k, l,) m, and the corners of the domain feature the vertex lies
// on: 1 for a vertex (side1), 2 for an edge (side2), 3 for a triangle (side3),
// none for a Voronoi vertex (side4). The vertex is where the feature meets the
// bisectors of i with each seed between i and m: as many equations as the
// feature has free directions.
#pragma once

#include <cmath>
#include <functional>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include <glm/glm.hpp>

#include "test_support.h"

#include "exact.h"
#include "kernel/predicate_counts.h"

struct Site
{
    glm::dvec3 p;
    double w = 0.0;
};

struct Case
{
    std::vector<Site> s;       // i first, m last, the bisector seeds between
    std::vector<glm::dvec3> q; // the feature's corners; empty for side4
};

// --- the construction, in any number type ---
//
// T is mpq_class for the oracle and double for the naive version. Every
// intermediate is stored in a T, never in `auto`: gmpxx's expression
// templates would otherwise keep references to temporaries.

template <class T> struct V3
{
    T x, y, z;
};

template <class T> V3<T> to_v3(const glm::dvec3 &v)
{
    return {T(v.x), T(v.y), T(v.z)};
}

template <class T> T dot(const V3<T> &a, const V3<T> &b)
{
    T r = a.x * b.x;
    r += a.y * b.y;
    r += a.z * b.z;
    return r;
}

template <class T> V3<T> minus(const V3<T> &a, const V3<T> &b)
{
    V3<T> r{a.x, a.y, a.z};
    r.x -= b.x;
    r.y -= b.y;
    r.z -= b.z;
    return r;
}

// pi(x) = |x - p|^2 - w, the power distance, by its definition.
template <class T> T power(const V3<T> &x, const Site &s)
{
    const V3<T> d = minus(x, to_v3<T>(s.p));
    T r = dot(d, d);
    r -= T(s.w);
    return r;
}

inline int sign_of(const mpq_class &v)
{
    return sgn(v);
}
inline int sign_of(double v)
{
    return (v > 0) - (v < 0);
}

// Roughly how large a value is, to choose pivots by.
inline double magnitude(double v)
{
    return std::abs(v);
}
inline double magnitude(const mpq_class &v)
{
    return std::abs(v.get_d());
}

// The vertex of a Case, or nothing if its system is singular (the feature is
// parallel to a bisector, or the seeds are coplanar): then the vertex does
// not exist, and a predicate must not be asked about it.
//
// Unknowns: the coordinates of x along the feature, x = q0 + sum t_k (q_k -
// q0); or x itself, for side4. Equations: pi_i(x) = pi_a(x) for each
// bisector seed a, which is linear in x because the |x|^2 cancels:
//
//   2 x . (p_a - p_i) = |p_a|^2 - |p_i|^2 + w_i - w_a
//
// Solved by Gaussian elimination, exact in rationals.
template <class T> std::optional<V3<T>> construct(const Case &c)
{
    const Site &i = c.s.front();
    const std::size_t n = c.s.size() - 2; // equations = unknowns
    V3<T> origin{T(0), T(0), T(0)};
    std::vector<V3<T>> dirs;
    if (c.q.empty())
    {
        dirs = {{T(1), T(0), T(0)}, {T(0), T(1), T(0)}, {T(0), T(0), T(1)}};
    }
    else
    {
        origin = to_v3<T>(c.q[0]);
        for (std::size_t k = 1; k < c.q.size(); ++k)
            dirs.push_back(minus(to_v3<T>(c.q[k]), origin));
    }
    if (dirs.size() != n)
        return std::nullopt; // a malformed Case; the drivers never make one

    // Rows A t = b, one per bisector seed.
    std::vector<std::vector<T>> A(n, std::vector<T>(n));
    std::vector<T> b(n);
    const V3<T> pi = to_v3<T>(i.p);
    for (std::size_t r = 0; r < n; ++r)
    {
        const Site &a = c.s[r + 1];
        const V3<T> pa = to_v3<T>(a.p);
        V3<T> normal = minus(pa, pi);
        normal.x *= T(2);
        normal.y *= T(2);
        normal.z *= T(2);
        T rhs = dot(pa, pa);
        rhs -= dot(pi, pi);
        rhs += T(i.w);
        rhs -= T(a.w);
        rhs -= dot(normal, origin);
        for (std::size_t k = 0; k < n; ++k)
            A[r][k] = dot(normal, dirs[k]);
        b[r] = rhs;
    }
    // Elimination with the largest pivot (for doubles; any nonzero one is
    // exact for rationals).
    for (std::size_t col = 0; col < n; ++col)
    {
        std::size_t piv = col;
        for (std::size_t r = col + 1; r < n; ++r)
            if (magnitude(A[r][col]) > magnitude(A[piv][col]))
                piv = r;
        if (A[piv][col] == 0)
            return std::nullopt;
        std::swap(A[piv], A[col]);
        std::swap(b[piv], b[col]);
        for (std::size_t r = col + 1; r < n; ++r)
        {
            const T f = A[r][col] / A[col][col];
            for (std::size_t k = col; k < n; ++k)
                A[r][k] -= f * A[col][k];
            b[r] -= f * b[col];
        }
    }
    std::vector<T> t(n);
    for (std::size_t r = n; r-- > 0;)
    {
        T v = b[r];
        for (std::size_t k = r + 1; k < n; ++k)
            v -= A[r][k] * t[k];
        t[r] = v / A[r][r];
    }
    V3<T> x = origin;
    for (std::size_t k = 0; k < n; ++k)
    {
        x.x += t[k] * dirs[k].x;
        x.y += t[k] * dirs[k].y;
        x.z += t[k] * dirs[k].z;
    }
    return x;
}

// The sign of pi_i(x) - pi_m(x) at the Case's vertex: negative on i's side.
// Nothing if the vertex does not exist.
template <class T> std::optional<int> construct_and_compare(const Case &c)
{
    const std::optional<V3<T>> x = construct<T>(c);
    if (!x)
        return std::nullopt;
    T d = power(*x, c.s.front());
    d -= power(*x, c.s.back());
    return sign_of(d);
}

inline std::optional<int> oracle(const Case &c)
{
    return construct_and_compare<mpq_class>(c);
}

// The floating-point way: build the vertex in doubles, then compare. Only for
// showing that an input set is hard; never an answer.
inline std::optional<int> naive(const Case &c)
{
    return construct_and_compare<double>(c);
}

// --- inputs ---

// Coordinates with random sign, full significand and exponent in [-20, 20]:
// cancellation and absorption at every scale, and products of six of them
// stay far from overflow and underflow (the filter's and DynamicExpansion's
// ranges).
inline double coordinate(std::mt19937_64 &rng)
{
    return random_double(rng, -20, 20);
}

inline Site random_site(std::mt19937_64 &rng)
{
    Site s{{coordinate(rng), coordinate(rng), coordinate(rng)}, 0.0};
    if (rng() & 1) // half of them weighted, half plain Voronoi
        s.w = coordinate(rng);
    return s;
}

inline glm::dvec3 random_point(std::mt19937_64 &rng)
{
    return {coordinate(rng), coordinate(rng), coordinate(rng)};
}

// Small integers: coordinates in [-2, 2], weights in [-2, 2] or zero. On so
// coarse a lattice, ties are common -- vertices exactly on a bisector, five
// seeds on one sphere -- and every input is exact in doubles.
inline Site lattice_site(std::mt19937_64 &rng)
{
    auto c = [&]() { return double(int(rng() % 5) - 2); };
    Site s{{c(), c(), c()}, 0.0};
    if (rng() % 3 == 0)
        s.w = c();
    return s;
}

inline glm::dvec3 lattice_point(std::mt19937_64 &rng)
{
    auto c = [&]() { return double(int(rng() % 5) - 2); };
    return {c(), c(), c()};
}

inline Case make_case(std::mt19937_64 &rng, std::size_t sites,
                      std::size_t corners, bool lattice)
{
    Case c;
    for (std::size_t k = 0; k < sites; ++k)
        c.s.push_back(lattice ? lattice_site(rng) : random_site(rng));
    for (std::size_t k = 0; k < corners; ++k)
        c.q.push_back(lattice ? lattice_point(rng) : random_point(rng));
    return c;
}

// The same Case moved by an exact offset: every coordinate plus `d`, weights
// unchanged. Differences of power distances do not change under a common
// translation, so neither do the answers -- but the cancellation the
// predicate must survive grows with d.
inline Case translated(Case c, double d)
{
    for (Site &s : c.s)
        s.p += glm::dvec3(d);
    for (glm::dvec3 &q : c.q)
        q += glm::dvec3(d);
    return c;
}

// One nonzero coordinate of one input moved by one ulp, up or down: from an
// exact tie, a vertex a hair off its bisector, where doubles lose the sign.
//
// Never a zero coordinate: one ulp up from 0 is the smallest subnormal,
// 2^-1074, whose square underflows to zero. That is outside the range where
// the error-free transformations are exact (eft.h: products above 2^-969),
// so no predicate built on them can be exact there either. Nothing is nudged
// if every coordinate is zero.
inline Case nudged(Case c, std::mt19937_64 &rng)
{
    std::vector<double *> nonzero;
    for (Site &s : c.s)
        for (int k = 0; k < 3; ++k)
            if (s.p[k] != 0)
                nonzero.push_back(&s.p[k]);
    for (glm::dvec3 &q : c.q)
        for (int k = 0; k < 3; ++k)
            if (q[k] != 0)
                nonzero.push_back(&q[k]);
    if (nonzero.empty())
        return c;
    double &x = *nonzero[rng() % nonzero.size()];
    x = std::nextafter(x, (rng() & 1) ? INFINITY : -INFINITY);
    return c;
}

inline std::string describe(const Case &c)
{
    auto v = [](const glm::dvec3 &p)
    { return "(" + show(p.x) + ", " + show(p.y) + ", " + show(p.z) + ")"; };
    std::string r = "seeds";
    for (const Site &s : c.s)
        r += " " + v(s.p) + " w " + show(s.w);
    if (!c.q.empty())
    {
        r += "; corners";
        for (const glm::dvec3 &q : c.q)
            r += " " + v(q);
    }
    return r;
}

struct SideTally
{
    long failures = 0;
    std::string first;

    void add(bool ok, const Case &c)
    {
        if (!ok && failures++ == 0)
            first = describe(c);
    }
    void report(const std::string &what) const
    {
        check(failures == 0, what.c_str(), failures ? first.c_str() : nullptr);
    }
};

// --- the driver: what every side predicate must pass ---
//
// `call` evaluates the predicate on a Case, `counts` is its counter. Five
// families of input, each against the oracle, skipping vertices that do not
// exist:
//
//   random     generic positions over many scales; the filter should decide
//              nearly all of them
//   lattice    small integers, many exact ties: the answer must be 0 there
//              (no perturbation yet -- that is step 7)
//   far        the lattice moved by 2^20 and by 10^6 + 0.5: the same answers,
//              ties included, under heavy cancellation
//   nudged     a tie moved by one ulp: tiny nonzero answers; these must reach
//              the exact path, and doubles must get some of them wrong
//
// The thresholds were measured on a reference implementation when the suites
// were written (2026-10-01), and set below what it reached:
//
//                filter on random   lattice ties   nudged reach exact
//   side1        100%               4.8%           89%
//   side2        100%               3.1%           92%
//   side3        98.1%              2.8%           91%
//   side4        99.98%             4.0%           93%
//
// Not every nudged tie needs the exact path: a nudge can move a value far
// enough from zero, relative to its error bound, for the filter to certify it.
// And constructing the vertex in doubles got 66% to 77% of them wrong.
struct DriverLimits
{
    double random_filtered; // at least this fraction of random cases
    double nudged_exact;    // at least this fraction of nudged cases
};

inline void check_predicate(const std::string &name, std::size_t sites,
                            std::size_t corners,
                            const std::function<int(const Case &)> &call,
                            PredicateCounts &counts, DriverLimits limits,
                            std::uint64_t seed)
{
    std::mt19937_64 rng(seed);

    // random
    {
        counts = PredicateCounts{};
        SideTally t;
        long n = 0;
        for (int it = 0; it < 50000; ++it)
        {
            const Case c = make_case(rng, sites, corners, false);
            const std::optional<int> want = oracle(c);
            if (!want)
                continue;
            ++n;
            t.add(call(c) == *want, c);
        }
        t.report(name + ": random inputs, the exact sign");
        check(counts.filtered + counts.exact == std::uint64_t(n),
              (name + ": every call is counted once").c_str());
        check(double(counts.filtered) >= limits.random_filtered * double(n),
              (name + ": the filter decides nearly all random inputs").c_str());
    }

    // lattice, and far
    {
        SideTally t, far20, far6;
        long ties = 0, n = 0;
        for (int it = 0; it < 50000; ++it)
        {
            const Case c = make_case(rng, sites, corners, true);
            const std::optional<int> want = oracle(c);
            if (!want)
                continue;
            ++n;
            ties += *want == 0;
            t.add(call(c) == *want, c);
            far20.add(call(translated(c, 0x1p20)) == *want, c);
            far6.add(call(translated(c, 1e6 + 0.5)) == *want, c);
        }
        t.report(name + ": lattice inputs, ties give 0");
        far20.report(name + ": the lattice moved by 2^20, the same answers");
        far6.report(name +
                    ": the lattice moved by 1e6 + 0.5, the same answers");
        check(ties >= n / 50, (name + ": the lattice has ties (at least 2%) "
                                      "-- or the test is too easy")
                                  .c_str());
    }

    // nudged ties
    {
        counts = PredicateCounts{};
        SideTally t;
        long n = 0, naive_wrong = 0;
        for (int it = 0; it < 200000 && n < 20000; ++it)
        {
            const Case tie = make_case(rng, sites, corners, true);
            const std::optional<int> at_tie = oracle(tie);
            if (!at_tie || *at_tie != 0)
                continue;
            const Case c = nudged(tie, rng);
            const std::optional<int> want = oracle(c);
            if (!want || *want == 0)
                continue;
            ++n;
            t.add(call(c) == *want, c);
            const std::optional<int> guess = naive(c);
            naive_wrong += !guess || *guess != *want;
        }
        t.report(name + ": ties nudged by one ulp, the exact sign");
        check(n >= 1000, (name + ": enough nudged ties to test").c_str());
        check(double(counts.exact) >= limits.nudged_exact * double(n),
              (name + ": nudged ties reach the exact path").c_str());
        check(naive_wrong > 0,
              (name + ": constructing the vertex in doubles gets some nudged "
                      "ties wrong (or the test is too easy)")
                  .c_str());
    }
}
