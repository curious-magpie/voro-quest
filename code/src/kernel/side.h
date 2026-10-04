#pragma once
#include "kernel/bounded.h"
#include "kernel/dynamic_expansion.h"
#include "kernel/predicate_counts.h"
#include <glm/ext/vector_double3.hpp>

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
