#pragma once
#include "kernel/bounded.h"
#include "kernel/dynamic_expansion.h"
#include "kernel/predicate_counts.h"
#include <glm/ext/vector_double3.hpp>
#include <utility>

inline PredicateCounts &side1_counts()
{
    thread_local PredicateCounts counts;
    return counts;
}

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

inline PredicateCounts &side2_counts()
{
    thread_local PredicateCounts counts;
    return counts;
}

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
