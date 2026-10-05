// The side predicates: which side of a bisector a vertex of a clipped piece
// is on, decided from the input numbers alone.
//
// Clipping a cell of seed i by the bisector of i and a new seed m keeps the
// part of the piece where i's power distance is the smaller one. Each vertex
// of the piece is asked one question: is pi_i(x) < pi_m(x) there? With the
// power distance pi_a(x) = |x - p_a|^2 - w_a, the answer is the sign of
//
//   f_im(x) = pi_i(x) - pi_m(x)
//
//   -1   x is on i's side: it stays
//   +1   x is on m's side: it is cut off
//    0   x is exactly on the bisector (a tie; step 7 will break it)
//
// The |x|^2 in the two power distances cancels, so every f_ia is affine in x.
// That one fact is what all four predicates are built on.
//
// --- the rule: no constructed coordinates ---
//
// A vertex x is never computed. It is known only by what defines it (design,
// section 5.5), and each predicate reads that definition:
//
//   side1   x is a domain vertex q                         degree 2
//   side2   x is where the edge q0 q1 crosses bisector ij   degree 4 / 2
//   side3   x is where the triangle q0 q1 q2 crosses
//           bisectors ij and ik                            degree 6 / 4
//   side4   x is the Voronoi vertex of i, j, k, l           degree 5 / 3
//
// Constructing x in doubles and then comparing the two distances is what a
// floating-point Voronoi code does, and it is wrong near a tie: on the
// tests' ties moved by one ulp, it got 65% to 76% of the signs wrong
// (tests/measure_side.cpp, 2026-10-05). Two
// cells that look at the same vertex through different roundings can then
// disagree, which is a crack in the diagram. Here each sign is the exact sign
// of a polynomial in the inputs, so every cell gets the same answer.
//
// side2 to side4 are ratios num / den. The sign of a ratio is the product of
// the two signs, so nothing is ever divided: each predicate returns
// sign(num) * sign(den). The denominator is zero exactly when the vertex's
// equations have no unique solution -- none, or infinitely many: an edge
// parallel to a bisector or lying in it, four coplanar seeds, two seeds at
// one position -- and a clipper never asks about such a vertex; the exact
// path asserts it.
//
// --- fast when easy, exact when not ---
//
// Each predicate is one formula, a template over the number type T, run
// twice at most, as orient3d is:
//
//   T = Bounded           the double value with a bound on its rounding
//                         error (bounded.h). If both signs are certain, or
//                         the numerator is certainly zero, that is the
//                         answer. On the tests' random inputs this decided
//                         98% to 100% of calls (tests/measure_side.cpp).
//   T = DynamicExpansion  the exact value (dynamic_expansion.h), when the
//                         filter cannot tell. Its length is set at runtime:
//                         a degree-6 formula would need millions of
//                         components as an Expansion<N> type.
//
// The counters (side1_counts() to side4_counts()) record which path decided
// each call, for the app's robustness window and for the tests.
//
// Exact within eft.h's range: products of input differences must not fall
// below 2^-969. A subnormal coordinate breaks that (its square underflows),
// so inputs are zero or of ordinary size.
//
// Not yet here: symbolic perturbation (design, section 6.4). A tie returns 0.
#pragma once
#include "kernel/bounded.h"
#include "kernel/dynamic_expansion.h"
#include "kernel/predicate_counts.h"
#include <cassert>
#include <glm/ext/vector_double3.hpp>
#include <utility>
#include <cstdint>

// How each side1 call was decided. One per predicate and per thread, like
// orient3d_counts(): thread-local, so parallel clipping needs no locks.
inline PredicateCounts &side1_counts()
{
    thread_local PredicateCounts counts;
    return counts;
}

// f_ia(q) = pi_i(q) - pi_a(q), in any number type T: negative when q is on
// i's side of the bisector of i and a. The building block of side1 to side3.
//
// Written in its translated form,
//
//   f_ia(q) = (p_a - p_i) . ((q - p_i) + (q - p_a)) + (w_a - w_i)
//
// which is the definition rearranged: with u = q - p_i and v = q - p_a,
// |u|^2 - |v|^2 = (u - v) . (u + v), and u - v = p_a - p_i. Every factor is
// then a difference of two inputs, so the filter's error bound scales with
// the distances involved, not with how far they are from the origin, and
// each difference is at most two components exactly.
//
// Each double is turned into a T before it is subtracted: T(pa.x) - T(pi.x),
// never T(pa.x - pi.x). The second rounds in double first, and Bounded would
// then call a rounded number exact, so the filter could certify a wrong sign.
template <class T>
T bisector_value(const glm::dvec3 &pi, double wi, const glm::dvec3 &pa,
                 double wa, const glm::dvec3 &q)
{
    const T dx = T(pa.x) - T(pi.x);
    const T dy = T(pa.y) - T(pi.y);
    const T dz = T(pa.z) - T(pi.z);

    const T sx = (T(q.x) - T(pi.x)) + (T(q.x) - T(pa.x));
    const T sy = (T(q.y) - T(pi.y)) + (T(q.y) - T(pa.y));
    const T sz = (T(q.z) - T(pi.z)) + (T(q.z) - T(pa.z));

    return dx * sx + dy * sy + dz * sz + (T(wa) - T(wi));
}

// Is the domain vertex q on i's side of the bisector of i and m? The sign of
// f_im(q): -1 on i's side, +1 on m's, 0 on the bisector.
inline int side1(const glm::dvec3 &pi, double wi, const glm::dvec3 &pm,
                 double wm, const glm::dvec3 &q)
{
    auto s = certain_sign(bisector_value<Bounded>(pi, wi, pm, wm, q));
    if (s)
    {
        ++side1_counts().filtered;
        return *s;
    }
    ++side1_counts().exact;
    return sign(bisector_value<DynamicExpansion>(pi, wi, pm, wm, q));
}

inline int side1_sos(std::uint32_t i, const glm::dvec3 &pi, double wi,
                     std::uint32_t m, const glm::dvec3 &pm, double wm,
                     const glm::dvec3 &q)
{
    const int s = side1(pi, wi, pm, wm, q);
    if (s != 0)
    {
        return s;
    }
    ++side1_counts().perturbed;
    assert(i != m);
    return i < m ? -1 : 1;
}

// --- side2: an edge crossing a bisector ---
inline PredicateCounts &side2_counts()
{
    thread_local PredicateCounts counts;
    return counts;
}

// The numerator and denominator of f_im at the point x where the edge q0 q1
// crosses the bisector of i and j. With m_r = f_im(q_r), j_r = f_ij(q_r):
//
//              m0 j1 - m1 j0
//   f_im(x) = ---------------
//                 j1 - j0
//
// Why: on the edge, x = q0 + t (q1 - q0), and an affine f is the same blend of
// its corner values, f(x) = (1 - t) f(q0) + t f(q1). The condition f_ij(x) =
// 0 fixes t = j0 / (j0 - j1), and putting that t into f_im gives the ratio.
// It is the 2x2 determinant [[m0, m1], [j0, j1]] over [[1, 1], [j0, j1]]
// (design, section 6.2): degree 4 over degree 2.
//
// Returned as a pair, so that one formula serves both number types; the
// caller unpacks it with `auto [num, den] = ...`.
template <class T>
std::pair<T, T> side2_terms(const glm::dvec3 &pi, double wi,
                            const glm::dvec3 &pj, double wj,
                            const glm::dvec3 &pm, double wm,
                            const glm::dvec3 &q0, const glm::dvec3 &q1)
{
    const T m0 = bisector_value<T>(pi, wi, pm, wm, q0); // f_im(q0)
    const T m1 = bisector_value<T>(pi, wi, pm, wm, q1); // f_im(q1)
    const T j0 = bisector_value<T>(pi, wi, pj, wj, q0); // f_ij(q0)
    const T j1 = bisector_value<T>(pi, wi, pj, wj, q1); // f_ij(q1)
    return {m0 * j1 - m1 * j0, j1 - j0};
}

// Is the point where the edge q0 q1 crosses bisector ij on i's side of the
// bisector of i and m? -1, +1 or 0, as side1. Precondition: the edge is not
// parallel to bisector ij (j0 != j1), so that the point exists.
//
// The filter is trusted in two cases: both signs certain, or a numerator that
// is certainly zero, which makes the answer 0 whatever the denominator (it is
// nonzero by the precondition). A certain numerator with an uncertain
// denominator is not enough: the sign of the product is then unknown.
inline int side2(const glm::dvec3 &pi, double wi, const glm::dvec3 &pj,
                 double wj, const glm::dvec3 &pm, double wm,
                 const glm::dvec3 &q0, const glm::dvec3 &q1)
{
    auto [num, den] = side2_terms<Bounded>(pi, wi, pj, wj, pm, wm, q0, q1);
    auto sn = certain_sign(num), sd = certain_sign(den);
    if (sn && *sn == 0)
    {
        ++side2_counts().filtered;
        return 0;
    }

    if (sn && sd)
    {
        ++side2_counts().filtered;
        return *sn * *sd;
    }

    ++side2_counts().exact;
    auto [en, ed] =
        side2_terms<DynamicExpansion>(pi, wi, pj, wj, pm, wm, q0, q1);
    assert(sign(ed) != 0);
    return sign(en) * sign(ed);
}

// --- side3: a triangle crossing two bisectors ---
inline PredicateCounts &side3_counts()
{
    thread_local PredicateCounts counts;
    return counts;
}

// The numerator and denominator of f_im at the point x where the triangle
// q0 q1 q2 crosses the bisectors of i with j and with k. With m_r, j_r, k_r
// the values of f_im, f_ij, f_ik at corner q_r:
//
//              | m0 m1 m2 |     | 1  1  1  |
//   f_im(x) =  | j0 j1 j2 |  /  | j0 j1 j2 |      degree 6 over degree 4
//              | k0 k1 k2 |     | k0 k1 k2 |
//
// Why: in barycentric coordinates x = l0 q0 + l1 q1 + l2 q2, with l0 + l1 +
// l2 = 1, and an affine f is f(x) = l0 f(q0) + l1 f(q1) + l2 f(q2). The point
// is fixed by three linear equations in the l's: they sum to 1, f_ij(x) = 0,
// f_ik(x) = 0, whose matrix is the denominator. Cramer's rule gives each l_r
// as the r-th cofactor of its first row over its determinant, so f_im(x) =
// sum m_r l_r is that determinant with its first row replaced by the m's.
//
// Both determinants are expanded along their first row, and below it they
// are the same: one set of three 2x2 minors c0, c1, c2 serves both. The
// middle cofactor has a minus sign, in both.
template <class T>
std::pair<T, T>
side3_terms(const glm::dvec3 &pi, double wi, const glm::dvec3 &pj, double wj,
            const glm::dvec3 &pk, double wk, const glm::dvec3 &pm, double wm,
            const glm::dvec3 &q0, const glm::dvec3 &q1, const glm::dvec3 &q2)
{
    const T m0 = bisector_value<T>(pi, wi, pm, wm, q0);
    const T m1 = bisector_value<T>(pi, wi, pm, wm, q1);
    const T m2 = bisector_value<T>(pi, wi, pm, wm, q2);

    const T j0 = bisector_value<T>(pi, wi, pj, wj, q0);
    const T j1 = bisector_value<T>(pi, wi, pj, wj, q1);
    const T j2 = bisector_value<T>(pi, wi, pj, wj, q2);

    const T k0 = bisector_value<T>(pi, wi, pk, wk, q0);
    const T k1 = bisector_value<T>(pi, wi, pk, wk, q1);
    const T k2 = bisector_value<T>(pi, wi, pk, wk, q2);

    const T c0 = j1 * k2 - j2 * k1;
    const T c1 = j0 * k2 - j2 * k0;
    const T c2 = j0 * k1 - j1 * k0;

    return {m0 * c0 - m1 * c1 + m2 * c2, c0 - c1 + c2};
}

// Is the point where the triangle q0 q1 q2 crosses bisectors ij and ik on
// i's side of the bisector of i and m? -1, +1 or 0. Precondition: the two
// bisectors meet the triangle's plane in one point (the denominator is not
// zero). The filter is trusted as in side2.
inline int side3(const glm::dvec3 &pi, double wi, const glm::dvec3 &pj,
                 double wj, const glm::dvec3 &pk, double wk,
                 const glm::dvec3 &pm, double wm, const glm::dvec3 &q0,
                 const glm::dvec3 &q1, const glm::dvec3 &q2)
{

    auto [num, den] =
        side3_terms<Bounded>(pi, wi, pj, wj, pk, wk, pm, wm, q0, q1, q2);
    auto sn = certain_sign(num), sd = certain_sign(den);
    if (sn && *sn == 0)
    {
        ++side3_counts().filtered;
        return 0;
    }

    if (sn && sd)
    {
        ++side3_counts().filtered;
        return *sn * *sd;
    }

    ++side3_counts().exact;
    auto [en, ed] = side3_terms<DynamicExpansion>(pi, wi, pj, wj, pk, wk, pm,
                                                  wm, q0, q1, q2);
    assert(sign(ed) != 0);
    return sign(en) * sign(ed);
}

// --- side4: a Voronoi vertex ---
inline PredicateCounts &side4_counts()
{
    thread_local PredicateCounts counts;
    return counts;
}

// One seed's row of side4's determinants: d = p_a - p_i, and c = |d|^2 + w_i
// - w_a. A class template, so Row<Bounded> and Row<DynamicExpansion> come
// from this one definition.
template <class T> struct Row
{
    T x, y, z, c;
};

// The 3x3 determinant of three rows' d parts (x, y, z; c is not read),
// expanded along the first row. T is deduced from the arguments.
template <class T> T det3(const Row<T> &a, const Row<T> &b, const Row<T> &c)
{
    return a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) +
           a.z * (b.x * c.y - b.y * c.x);
}

// The numerator and denominator of f_im at the Voronoi vertex x of i, j, k,
// l: the point with equal power distance to all four.
//
// In coordinates relative to seed i, y = x - p_i and d_a = p_a - p_i, the
// definition expands to f_ia(x) = 2 d_a . y - (|d_a|^2 + w_i - w_a). With
// z = 2y the 2 goes away:
//
//   f_ia(x) = d_a . z - c_a,    c_a = |d_a|^2 + w_i - w_a
//
// and x is where d_a . z = c_a for a = j, k, l. Note the weights: w_i - w_a
// here, the opposite order to bisector_value's, because c_a sits on the
// other side of the equation.
//
// Take det4, the 4x4 determinant with rows (d_a, c_a) for a = j, k, l, m.
// Subtracting z.x, z.y, z.z times the first three columns from the last does
// not change it, and turns the last column into c_a - d_a . z = -f_ia(x):
// zero for j, k, l, and -f_im(x) for m. Expanded along that column, only the
// bottom entry is left, so
//
//   det4 = -f_im(x) det3(d_j; d_k; d_l),   f_im(x) = -det4 / det3
//
// degree 5 over degree 3 (design, section 6.2 writes u_a = 2 d_a, which
// scales both determinants by 8 and changes no sign). det4 itself is
// computed by that same expansion along the c column, cofactor signs
// - + - + down it, each minor being det3 of the other three rows in their
// order; the last minor is the denominator, so it is computed once.
//
// The rows are built by a lambda that reads p_i and w_i from this function
// ([&]): one place for the translated d and c, used for four seeds.
template <class T>
std::pair<T, T>
side4_terms(const glm::dvec3 &pi, double wi, const glm::dvec3 &pj, double wj,
            const glm::dvec3 &pk, double wk, const glm::dvec3 &pl, double wl,
            const glm::dvec3 &pm, double wm)
{
    auto row = [&](const glm::dvec3 &pa, double wa)
    {
        const T dx = T(pa.x) - T(pi.x);
        const T dy = T(pa.y) - T(pi.y);
        const T dz = T(pa.z) - T(pi.z);
        return Row<T>{dx, dy, dz,
                      dx * dx + dy * dy + dz * dz + (T(wi) - T(wa))};
    };
    const Row<T> rj = row(pj, wj), rk = row(pk, wk), rl = row(pl, wl),
                 rm = row(pm, wm);
    const T den = det3(rj, rk, rl);
    const T det4 = rk.c * det3(rj, rl, rm) - rj.c * det3(rk, rl, rm) -
                   rl.c * det3(rj, rk, rm) + rm.c * den;
    return {-det4, den};
}

// Is the Voronoi vertex of i, j, k, l on i's side of the bisector of i and
// m? -1, +1, or 0 when the five seeds lie on one power sphere. Precondition:
// i, j, k, l are not coplanar (det3 is not zero), so that the vertex exists.
// The filter is trusted as in side2.
inline int side4(const glm::dvec3 &pi, double wi, const glm::dvec3 &pj,
                 double wj, const glm::dvec3 &pk, double wk,
                 const glm::dvec3 &pl, double wl, const glm::dvec3 &pm,
                 double wm)
{

    auto [num, den] =
        side4_terms<Bounded>(pi, wi, pj, wj, pk, wk, pl, wl, pm, wm);
    auto sn = certain_sign(num), sd = certain_sign(den);
    if (sn && *sn == 0)
    {
        ++side4_counts().filtered;
        return 0;
    }

    if (sn && sd)
    {
        ++side4_counts().filtered;
        return *sn * *sd;
    }

    ++side4_counts().exact;
    auto [en, ed] =
        side4_terms<DynamicExpansion>(pi, wi, pj, wj, pk, wk, pl, wl, pm, wm);
    assert(sign(ed) != 0);
    return sign(en) * sign(ed);
}
