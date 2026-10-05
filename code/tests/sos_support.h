// What the perturbed side-predicate suites share: the perturbed oracle, its
// own check, and the driver.
//
// Simulation of Simplicity (design, section 6.4): every seed's weight w_a is
// read as w_a + eps_a, where the eps are infinitesimals ordered by seed index,
// eps_0 >> eps_1 >> eps_2 >> ... -- the smaller the index, the infinitely
// larger its eps. Under that imaginary input no side predicate is ever zero,
// and all of them answer questions about the same input, so they agree.
//
// The oracle follows the definition, not the predicates' formulas. The vertex
// is solved exactly, as in side_support.h, but with every weight carrying its
// eps symbolically. The system's matrix holds no weights, so only its right
// side carries eps, and the solution, and then pi_i(x) - pi_m(x), are linear
// in the eps:
//
//   F = F_0 + sum over seeds a of F_a eps_a
//
// (the |x|^2 terms cancel exactly, eps^2 terms with them). Its sign is that
// of F_0 if F_0 is nonzero, and otherwise that of the first nonzero F_a,
// taking the seeds by increasing index.
//
// That oracle is itself checked, by doing what the eps stand for: give each
// seed a real, tiny extra weight, delta^(rank + 1) with delta = 2^-100 and
// rank its position in index order, and solve that input exactly. On the
// small lattice inputs where it is used, delta is far smaller than any ratio
// of the coefficients, so the two must agree.
#pragma once

#include <algorithm>
#include <numeric>

#include "side_support.h"

// Distinct random ids for the Case's seeds, drawn from [0, 64): the order of
// the arguments says nothing about the order of the indices.
inline Case with_ids(Case c, std::mt19937_64 &rng)
{
    std::vector<std::uint32_t> pool(64);
    std::iota(pool.begin(), pool.end(), 0u);
    std::shuffle(pool.begin(), pool.end(), rng);
    for (std::size_t k = 0; k < c.s.size(); ++k)
        c.s[k].id = pool[k];
    return c;
}

// --- the perturbed value, linear in the eps ---

// a constant plus one coefficient per seed of the Case, by position in c.s.
struct Linear
{
    mpq_class c;
    std::vector<mpq_class> e;
};

inline Linear linear_zero(std::size_t seeds)
{
    return Linear{mpq_class(0), std::vector<mpq_class>(seeds, mpq_class(0))};
}

// a -= f * b
inline void sub_scaled(Linear &a, const mpq_class &f, const Linear &b)
{
    mpq_class t = f * b.c;
    a.c -= t;
    for (std::size_t k = 0; k < a.e.size(); ++k)
    {
        t = f * b.e[k];
        a.e[k] -= t;
    }
}

// pi_i(x) - pi_m(x) at the Case's vertex, with every weight w_a + eps_a, as a
// Linear; nothing if the vertex does not exist. The unknowns and equations
// are those of construct() in side_support.h:
//
//   2 x . (p_a - p_i) = |p_a|^2 - |p_i|^2 + (w_i + eps_i) - (w_a + eps_a)
inline std::optional<Linear> perturbed_value(const Case &c)
{
    const std::size_t seeds = c.s.size(), n = seeds - 2;
    using Q = mpq_class;
    const V3<Q> pi = to_v3<Q>(c.s.front().p);
    V3<Q> origin{Q(0), Q(0), Q(0)};
    std::vector<V3<Q>> dirs;
    if (c.q.empty())
        dirs = {{Q(1), Q(0), Q(0)}, {Q(0), Q(1), Q(0)}, {Q(0), Q(0), Q(1)}};
    else
    {
        origin = to_v3<Q>(c.q[0]);
        for (std::size_t k = 1; k < c.q.size(); ++k)
            dirs.push_back(minus(to_v3<Q>(c.q[k]), origin));
    }

    std::vector<std::vector<Q>> A(n, std::vector<Q>(n));
    std::vector<Linear> b(n, linear_zero(seeds));
    for (std::size_t r = 0; r < n; ++r)
    {
        const Site &a = c.s[r + 1];
        const V3<Q> pa = to_v3<Q>(a.p);
        V3<Q> normal = minus(pa, pi);
        normal.x *= 2;
        normal.y *= 2;
        normal.z *= 2;
        Q rhs = dot(pa, pa);
        rhs -= dot(pi, pi);
        rhs += Q(c.s.front().w);
        rhs -= Q(a.w);
        rhs -= dot(normal, origin);
        b[r].c = rhs;
        b[r].e[0] += 1;     // + eps_i
        b[r].e[r + 1] -= 1; // - eps_a
        for (std::size_t k = 0; k < n; ++k)
            A[r][k] = dot(normal, dirs[k]);
    }
    for (std::size_t col = 0; col < n; ++col)
    {
        std::size_t piv = col;
        while (piv < n && A[piv][col] == 0)
            ++piv;
        if (piv == n)
            return std::nullopt;
        std::swap(A[piv], A[col]);
        std::swap(b[piv], b[col]);
        for (std::size_t r = col + 1; r < n; ++r)
        {
            const Q f = A[r][col] / A[col][col];
            for (std::size_t k = col; k < n; ++k)
            {
                const Q t = f * A[col][k];
                A[r][k] -= t;
            }
            sub_scaled(b[r], f, b[col]);
        }
    }
    std::vector<Linear> t(n, linear_zero(seeds));
    for (std::size_t r = n; r-- > 0;)
    {
        Linear v = b[r];
        for (std::size_t k = r + 1; k < n; ++k)
            sub_scaled(v, A[r][k], t[k]);
        const Q inv = 1 / A[r][r];
        v.c *= inv;
        for (Q &e : v.e)
            e *= inv;
        t[r] = v;
    }

    // x = origin + sum t_k dirs_k, then F = 2 x . (p_m - p_i) - (|p_m|^2 -
    // |p_i|^2 + w_i - w_m) - eps_i + eps_m.
    const Site &m = c.s.back();
    const V3<Q> pm = to_v3<Q>(m.p);
    V3<Q> g = minus(pm, pi);
    g.x *= 2;
    g.y *= 2;
    g.z *= 2;
    Linear F = linear_zero(seeds);
    F.c = dot(g, origin);
    for (std::size_t k = 0; k < n; ++k)
    {
        const Q gk = dot(g, dirs[k]); // F gains gk * t_k
        sub_scaled(F, -gk, t[k]);
    }
    Q k0 = dot(pm, pm);
    k0 -= dot(pi, pi);
    k0 += Q(c.s.front().w);
    k0 -= Q(m.w);
    F.c -= k0;
    F.e[0] -= 1;
    F.e[seeds - 1] += 1;
    return F;
}

// The sign under the perturbation: F_0's, or the first nonzero coefficient's
// in increasing seed index. 0 only if F is identically zero, which the
// suites check never happens for distinct seeds.
inline int perturbed_sign(const Linear &F, const Case &c)
{
    if (F.c != 0)
        return sgn(F.c);
    std::vector<std::size_t> order(c.s.size());
    std::iota(order.begin(), order.end(), std::size_t(0));
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
              { return c.s[a].id < c.s[b].id; });
    for (std::size_t k : order)
        if (F.e[k] != 0)
            return sgn(F.e[k]);
    return 0;
}

inline std::optional<int> perturbed_oracle(const Case &c)
{
    const std::optional<Linear> F = perturbed_value(c);
    if (!F)
        return std::nullopt;
    return perturbed_sign(*F, c);
}

// --- the oracle's own check: a real, tiny perturbation ---
//
// pi_i(x) - pi_m(x) with weights w_a + delta^(rank_a + 1), solved exactly by
// the same elimination with rational weights, and compared from the power
// distances' definition.
inline std::optional<int> concretely_perturbed(const Case &c)
{
    using Q = mpq_class;
    const std::size_t seeds = c.s.size(), n = seeds - 2;
    std::vector<std::size_t> order(seeds);
    std::iota(order.begin(), order.end(), std::size_t(0));
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
              { return c.s[a].id < c.s[b].id; });
    std::vector<Q> w(seeds);
    const Q delta(mpz_class(1), mpz_class(1) << 100);
    Q eps = delta;
    for (std::size_t k : order)
    {
        w[k] = Q(c.s[k].w);
        w[k] += eps;
        eps *= delta;
    }

    const V3<Q> pi = to_v3<Q>(c.s.front().p);
    V3<Q> origin{Q(0), Q(0), Q(0)};
    std::vector<V3<Q>> dirs;
    if (c.q.empty())
        dirs = {{Q(1), Q(0), Q(0)}, {Q(0), Q(1), Q(0)}, {Q(0), Q(0), Q(1)}};
    else
    {
        origin = to_v3<Q>(c.q[0]);
        for (std::size_t k = 1; k < c.q.size(); ++k)
            dirs.push_back(minus(to_v3<Q>(c.q[k]), origin));
    }
    std::vector<std::vector<Q>> A(n, std::vector<Q>(n));
    // Filled with explicit zeros: `b(n)` alone draws a -Wfree-nonheap-object
    // false positive from GCC 16 when n can be 0 (side1).
    std::vector<Q> b(n, Q(0));
    for (std::size_t r = 0; r < n; ++r)
    {
        const V3<Q> pa = to_v3<Q>(c.s[r + 1].p);
        V3<Q> normal = minus(pa, pi);
        normal.x *= 2;
        normal.y *= 2;
        normal.z *= 2;
        Q rhs = dot(pa, pa);
        rhs -= dot(pi, pi);
        rhs += w[0];
        rhs -= w[r + 1];
        rhs -= dot(normal, origin);
        b[r] = rhs;
        for (std::size_t k = 0; k < n; ++k)
            A[r][k] = dot(normal, dirs[k]);
    }
    for (std::size_t col = 0; col < n; ++col)
    {
        std::size_t piv = col;
        while (piv < n && A[piv][col] == 0)
            ++piv;
        if (piv == n)
            return std::nullopt;
        std::swap(A[piv], A[col]);
        std::swap(b[piv], b[col]);
        for (std::size_t r = col + 1; r < n; ++r)
        {
            const Q f = A[r][col] / A[col][col];
            for (std::size_t k = col; k < n; ++k)
            {
                const Q t = f * A[col][k];
                A[r][k] -= t;
            }
            const Q t = f * b[col];
            b[r] -= t;
        }
    }
    std::vector<Q> t(n);
    for (std::size_t r = n; r-- > 0;)
    {
        Q v = b[r];
        for (std::size_t k = r + 1; k < n; ++k)
        {
            const Q u = A[r][k] * t[k];
            v -= u;
        }
        t[r] = v / A[r][r];
    }
    V3<Q> x = origin;
    for (std::size_t k = 0; k < n; ++k)
    {
        Q u = t[k] * dirs[k].x;
        x.x += u;
        u = t[k] * dirs[k].y;
        x.y += u;
        u = t[k] * dirs[k].z;
        x.z += u;
    }
    const V3<Q> di = minus(x, pi), dm = minus(x, to_v3<Q>(c.s.back().p));
    Q d = dot(di, di);
    d -= w[0];
    d -= dot(dm, dm);
    d += w[seeds - 1];
    return sgn(d);
}

// The Case scaled by s: coordinates times s, weights times s^2. Every power
// distance scales by s^2, so ties stay ties and every perturbed answer stays
// the same. With s of 24 significant bits, s^2 is still an exact double and
// so is every scaled lattice input -- but products of three or more of them
// no longer fit in 53 bits. On the plain lattice every intermediate value is
// a small integer, and doubles would be exact by luck; here they are not.
inline Case scaled(Case c, double s)
{
    for (Site &st : c.s)
    {
        st.p *= s;
        st.w *= s * s;
    }
    for (glm::dvec3 &q : c.q)
        q *= s;
    return c;
}

// --- the driver: what every perturbed side predicate must pass ---
//
// `exact` is the unperturbed predicate (0 on a tie), `sos` the perturbed one,
// `counts` their shared counter: the perturbed one is expected to call the
// exact one and count a resolved tie in `perturbed`. On the lattice, where
// ties are common, and on random inputs:
//
//   oracle    the perturbed oracle agrees with a real tiny perturbation, and
//             with the unperturbed oracle wherever that is nonzero
//   never 0   sos never returns 0
//   agrees    sos equals exact wherever exact is nonzero
//   ties      at the ties, sos equals the perturbed oracle -- also with the
//             lattice moved by 2^20, scaled by a 24-bit factor (where doubles
//             round), unweighted and scaled by a 40-bit factor, and with the
//             ids moved in the same order
//   counted   counts.perturbed is exactly the number of ties
inline void check_perturbed(const std::string &name, std::size_t sites,
                            std::size_t corners,
                            const std::function<int(const Case &)> &exact,
                            const std::function<int(const Case &)> &sos,
                            PredicateCounts &counts, std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
    SideTally self, self0, never0, agrees, ties, far, scaled_ok, unweighted,
        relabelled;
    long n = 0, n_ties = 0, n_ties_u = 0;
    counts = PredicateCounts{};
    for (int it = 0; it < 30000; ++it)
    {
        const Case c = with_ids(make_case(rng, sites, corners, true), rng);
        const std::optional<int> want0 = oracle(c);
        if (!want0)
            continue;
        const std::optional<int> want = perturbed_oracle(c);
        ++n;
        n_ties += *want0 == 0;
        self.add(want && *want == concretely_perturbed(c), c);
        self0.add(want && *want != 0 && (*want0 == 0 || *want == *want0), c);

        const int e = exact(c), s = sos(c);
        never0.add(s != 0, c);
        agrees.add(e == 0 || s == e, c);
        ties.add(s == *want, c);
        far.add(sos(translated(c, 0x1p20)) == *want, c);
        // s = 0x1.aaaaaap+0: 24 significant bits, its last one set.
        scaled_ok.add(sos(scaled(c, 0x1.aaaaaap+0)) == *want, c);
        // Without weights, a wider factor keeps every input exact: 40 bits,
        // so even products of two inputs round in doubles. A different input
        // (the weights are gone), so its own oracle and its own ties.
        Case u = c;
        for (Site &st : u.s)
            st.w = 0;
        const std::optional<int> want_u = perturbed_oracle(u);
        if (want_u && oracle(u))
        {
            unweighted.add(sos(scaled(u, 0x1.aaaaaaaaaap+0)) == *want_u, u);
            n_ties_u += *oracle(u) == 0;
        }
        // Ids moved but kept in the same order: the same answers. (A
        // different order is a different perturbation; the ties check above
        // already draws a new order for every input.)
        Case shifted = c;
        for (Site &st : shifted.s)
            st.id = 3 * st.id + 7;
        relabelled.add(sos(shifted) == s, c);
    }
    self.report(name + ": the perturbed oracle agrees with a real tiny "
                       "perturbation");
    self0.report(name + ": the perturbed oracle is never 0, and keeps every "
                        "nonzero answer");
    never0.report(name + ": the perturbed predicate never returns 0");
    agrees.report(name + ": it equals the exact predicate wherever that is "
                         "not 0");
    ties.report(name + ": at the ties, it equals the perturbed oracle");
    far.report(name + ": the lattice moved by 2^20, the same answers");
    scaled_ok.report(name + ": the lattice scaled by a 24-bit factor, the "
                            "same answers");
    unweighted.report(name + ": unweighted, scaled by a 40-bit factor, the "
                             "perturbed oracle's answers");
    relabelled.report(name + ": ids moved in the same order, the same "
                             "answers");
    check(n_ties >= n / 50,
          (name + ": the lattice has ties to break (at least 2%)").c_str());
    // Each input was evaluated five times: exact, sos, far, scaled,
    // relabelled; the last four have the input's ties. The unweighted copy
    // has its own.
    check(counts.perturbed == std::uint64_t(4 * n_ties + n_ties_u),
          (name + ": every resolved tie is counted, and nothing else").c_str());

    // Random inputs have no ties: the perturbation never changes anything.
    SideTally random;
    for (int it = 0; it < 20000; ++it)
    {
        const Case c = with_ids(make_case(rng, sites, corners, false), rng);
        if (!oracle(c))
            continue;
        const int e = exact(c);
        random.add(e != 0 && sos(c) == e, c);
    }
    random.report(name + ": random inputs, no ties, the exact answer");
}
