// orient3d: which side of a plane a point is on -- orient2d, one dimension up.
//
// orient3d(a, b, c, d) is the sign of six times the signed volume of the tet
// a, b, c, d:
//
//   positive  d is on the side of the plane through a, b, c that a, b, c wind
//             counterclockwise around: the tet is positively oriented
//   negative  d is on the other side
//   zero      the four points are coplanar
//
// --- the convention, and why it is not Shewchuk's ---
//
// This is the orientation the rest of the code already uses: mesh-dev's
// TetMesh puts every tet in this order at build time, TetMesh::volume is
// positive for it, and the viewer's MeshView asks for tets "in positive order"
// in the same sense. Shewchuk's orient3d has the opposite sign -- positive
// when d is *below* the plane, "below" meaning the side a, b, c look clockwise
// from -- which would make every positive tet negative here. So the sign
// follows the code RVD runs with, not the paper, and
// tests/test_orient3d.cpp checks it against TetMesh itself: every tet of
// mesh-dev's Kuhn cube must come out positive.
//
// --- the formula ---
//
// Translated so that d is the origin, as orient2d translates by c and for the
// same reason: the filter's rounding error then scales with the size of the
// tet, not with how far it is from the origin (see orient2d.h). With
// A = a - d, B = b - d, C = c - d, the determinant is the triple product
//
//   A . (C x B) =   Ax (Cy Bz - Cz By)
//                 + Ay (Cz Bx - Cx Bz)
//                 + Az (Cx By - Cy Bx)
//
// C x B, not B x C: swapping the two is exactly the sign flip between this
// convention and Shewchuk's. On the unit tet -- a at the origin, b, c, d on
// the x, y, z axes -- it gives +1.
//
// --- one formula, three number types ---
//
// As in orient2d.h: orient3d_det is written once, and evaluated with Bounded
// for the filter, with Expansion<1> for the exact path, and with double where
// a rounded value is wanted. The exact capacities, worked out by the compiler,
// are
//
//   differences          Expansion<2>
//   Cy * Bz              Expansion<8>
//   Cy * Bz - Cz * By    Expansion<16>
//   Ax * (...)           Expansion<64>
//   the three terms      Expansion<192>
//
// 192 doubles, 1.5 KB on the stack -- the same bound as Shewchuk's exact
// orient3d. On random points the filter decides more than 99.9% of calls, so
// the exact path is reached by points a few ulps from coplanar, and by
// coplanar ones, nearly all of which need it (a rounded zero cannot be
// certified; only an exact zero difference, d on a, b or c, lets the filter
// certify one).
//
// --- range ---
//
// Exact while no product overflows or underflows (kernel/eft.h), and here
// three differences are multiplied, so the range is narrower than orient2d's:
// coordinates zero, or of magnitude between about 2^-250 and 2^330. Mesh and
// seed coordinates are far inside it, and the domain and seed set will check.
#pragma once

#include "kernel/bounded.h"
#include "kernel/expansion.h"
#include "kernel/predicate_counts.h"
#include <type_traits>

inline PredicateCounts &orient3d_counts()
{
    thread_local PredicateCounts counts;
    return counts;
}

// Six times the signed volume of the tet a, b, c, d, for any number type with
// -, * and +. The differences are named once and reused, so each is computed
// once; their type is whatever T's subtraction gives (Expansion<2> for
// Expansion<1> inputs).
template <class T>
auto orient3d_det(const T &ax, const T &ay, const T &az, const T &bx,
                  const T &by, const T &bz, const T &cx, const T &cy,
                  const T &cz, const T &dx, const T &dy, const T &dz)
{
    const auto Ax = ax - dx, Ay = ay - dy, Az = az - dz; // a - d
    const auto Bx = bx - dx, By = by - dy, Bz = bz - dz; // b - d
    const auto Cx = cx - dx, Cy = cy - dy, Cz = cz - dz; // c - d
    return Ax * (Cy * Bz - Cz * By) + Ay * (Cz * Bx - Cx * Bz) +
           Az * (Cx * By - Cy * Bx);
}

// The exact sign: +1 for a positively oriented tet, -1 for a negative one, 0
// for coplanar points. The filter first -- the same template, with Bounded --
// and the exact path only when it cannot decide; each call counts which of the
// two decided it, in orient3d_counts().
inline int orient3d(double ax, double ay, double az, double bx, double by,
                    double bz, double cx, double cy, double cz, double dx,
                    double dy, double dz)
{
    auto s = certain_sign(orient3d_det(Bounded(ax), Bounded(ay), Bounded(az),
                                       Bounded(bx), Bounded(by), Bounded(bz),
                                       Bounded(cx), Bounded(cy), Bounded(cz),
                                       Bounded(dx), Bounded(dy), Bounded(dz)));
    if (s)
    {
        ++orient3d_counts().filtered;
        return *s;
    }
    ++orient3d_counts().exact;
    using E = Expansion<1>;
    const auto det = orient3d_det(E(ax), E(ay), E(az), E(bx), E(by), E(bz),
                                  E(cx), E(cy), E(cz), E(dx), E(dy), E(dz));
    static_assert(std::is_same_v<decltype(det), const Expansion<192>>);
    return sign(det);
}
