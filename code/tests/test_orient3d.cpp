// orient3d (kernel/orient3d.h): which side of a plane a point is on, against
// exact rationals, and against mesh-dev's tet orientation.
//
// orient3d(a, b, c, d) is the sign of six times the signed volume of the tet
// a, b, c, d: positive when d is on the side of the plane abc that a, b, c wind
// counterclockwise around -- the convention of TetMesh and the viewer, and the
// opposite of Shewchuk's. The oracle is written independently, as the triple
// product (b - a) . ((c - a) x (d - a)) in rationals.
//
// As for orient2d: random points, points a few ulps off a plane where doubles
// fail (checked), points exactly on a plane where the answer is zero, the
// symmetries of a determinant, and the filter's counters. And one more: every
// tet of mesh-dev's Kuhn cube, which TetMesh::build orients positively, must
// come out positive here.
#include <cmath>
#include <random>
#include <string>

#include "test_support.h"

#include "exact.h"
#include "kernel/orient3d.h"
#include "mesh/tet_mesh.h"

namespace
{

struct P
{
    double x, y, z;
};

int oracle(const P &a, const P &b, const P &c, const P &d)
{
    const mpq_class bx = exact(b.x) - exact(a.x), by = exact(b.y) - exact(a.y),
                    bz = exact(b.z) - exact(a.z);
    const mpq_class cx = exact(c.x) - exact(a.x), cy = exact(c.y) - exact(a.y),
                    cz = exact(c.z) - exact(a.z);
    const mpq_class dx = exact(d.x) - exact(a.x), dy = exact(d.y) - exact(a.y),
                    dz = exact(d.z) - exact(a.z);
    const mpq_class v =
        bx * (cy * dz - cz * dy) + by * (cz * dx - cx * dz) + bz * (cx * dy - cy * dx);
    return sgn(v);
}

int exact_sign(const P &a, const P &b, const P &c, const P &d)
{
    return orient3d(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z);
}

int double_sign(const P &a, const P &b, const P &c, const P &d)
{
    const double v = orient3d_det(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z,
                                  d.x, d.y, d.z);
    return (v > 0) - (v < 0);
}

std::string describe(const P &a, const P &b, const P &c, const P &d)
{
    std::string s;
    for (const P *p : {&a, &b, &c, &d})
        s += "(" + show(p->x) + ", " + show(p->y) + ", " + show(p->z) + ") ";
    return s;
}

struct Tally
{
    long failures = 0;
    std::string first;

    void add(bool ok, const P &a, const P &b, const P &c, const P &d)
    {
        if (!ok && failures++ == 0)
            first = describe(a, b, c, d);
    }
    void report(const char *what) const
    {
        check(failures == 0, what, failures ? first.c_str() : nullptr);
    }
};

void test_random()
{
    orient3d_counts() = PredicateCounts{};
    Tally t;
    std::mt19937_64 rng(70);
    const long n = 100000;
    for (long i = 0; i < n; ++i)
    {
        P p[4];
        for (P &q : p)
            q = {random_double(rng, -30, 30), random_double(rng, -30, 30),
                 random_double(rng, -30, 30)};
        t.add(exact_sign(p[0], p[1], p[2], p[3]) == oracle(p[0], p[1], p[2], p[3]),
              p[0], p[1], p[2], p[3]);
    }
    t.report("orient3d: random points, the exact sign");
    const PredicateCounts c = orient3d_counts();
    check(c.filtered + c.exact == uint64_t(n),
          "orient3d: every call is counted once");
    check(c.filtered >= uint64_t(0.999 * n),
          "orient3d: the filter decides at least 99.9% of random points");
}

// The 3D counterpart of Kettner et al.'s grid: points one ulp apart around
// (0.5, 0.25, 0.5), against the plane z = x through three points far away.
void test_near_plane()
{
    orient3d_counts() = PredicateCounts{};
    Tally t;
    long double_wrong = 0;
    const P a{12, 0, 12}, b{24, 7, 24}, c{17, -5, 17};
    const double ulp = std::ldexp(1.0, -53);
    for (int i = 0; i < 128; ++i)
        for (int j = 0; j < 128; ++j)
        {
            const P d{0.5 + i * ulp, 0.25, 0.5 + j * ulp};
            const int want = oracle(a, b, c, d);
            t.add(exact_sign(a, b, c, d) == want, a, b, c, d);
            double_wrong += double_sign(a, b, c, d) != want;
        }
    t.report("orient3d: a grid one ulp apart near a plane, the exact sign");
    check(double_wrong > 0,
          "the double formula gets some of that grid wrong (or the test is "
          "too easy)");
    check(orient3d_counts().exact > 0, "the near-plane grid reaches the exact path");
}

// Exactly on a plane: d = a + s (b - a) + t (c - a) with integers.
void test_on_plane()
{
    orient3d_counts() = PredicateCounts{};
    Tally t;
    std::mt19937_64 rng(71);
    std::uniform_int_distribution<long> coord(-(1L << 30), 1L << 30);
    std::uniform_int_distribution<int> k(-20, 20);
    const long n = 20000;
    for (long i = 0; i < n; ++i)
    {
        const P a{double(coord(rng)), double(coord(rng)), double(coord(rng))};
        const P b{double(coord(rng)), double(coord(rng)), double(coord(rng))};
        const P c{double(coord(rng)), double(coord(rng)), double(coord(rng))};
        const int s = k(rng), u = k(rng);
        const P d{a.x + s * (b.x - a.x) + u * (c.x - a.x),
                  a.y + s * (b.y - a.y) + u * (c.y - a.y),
                  a.z + s * (b.z - a.z) + u * (c.z - a.z)};
        t.add(exact_sign(a, b, c, d) == 0, a, b, c, d);
    }
    t.report("orient3d: exactly coplanar points give zero");
    // Nearly all of them must go exact: a rounded zero cannot be certified.
    // Not all: when d lands on a, b or c, a difference is an exact zero, its
    // error is zero too, and the filter rightly certifies the zero itself.
    check(orient3d_counts().exact >= uint64_t(0.99 * n),
          "orient3d: coplanar points nearly always reach the exact path");
}

// A determinant's symmetries: an odd permutation of the points flips the sign,
// an even one keeps it. Far from the origin, where cancellation is worst.
void test_symmetry()
{
    Tally odd, even, far;
    std::mt19937_64 rng(72);
    for (int i = 0; i < 50000; ++i)
    {
        P p[4];
        for (P &q : p)
            q = {1e6 + random_double(rng, -20, -10), 1e6 + random_double(rng, -20, -10),
                 1e6 + random_double(rng, -20, -10)};
        const int s = exact_sign(p[0], p[1], p[2], p[3]);
        odd.add(exact_sign(p[1], p[0], p[2], p[3]) == -s &&
                    exact_sign(p[0], p[1], p[3], p[2]) == -s,
                p[0], p[1], p[2], p[3]);
        even.add(exact_sign(p[1], p[2], p[0], p[3]) == s &&
                     exact_sign(p[1], p[0], p[3], p[2]) == s,
                 p[0], p[1], p[2], p[3]);
        far.add(s == oracle(p[0], p[1], p[2], p[3]), p[0], p[1], p[2], p[3]);
    }
    odd.report("orient3d: an odd permutation flips the sign");
    even.report("orient3d: an even permutation keeps the sign");
    far.report("orient3d: points near 1e6, the exact sign");
}

// mesh-dev's convention, from mesh-dev's own code: TetMesh::build puts every
// tet of the Kuhn cube in positive order, so orient3d must call each positive.
void test_tetmesh_convention()
{
    TetMesh mesh;
    check(build_tet_cube(mesh, 4, 1.0), "mesh-dev's tet cube builds");
    long wrong = 0;
    for (uint32_t t = 0; t < mesh.tet_count(); ++t)
    {
        const uint32_t *v = mesh.tet(t);
        const glm::dvec3 &a = mesh.p(v[0]), &b = mesh.p(v[1]), &c = mesh.p(v[2]),
                         &d = mesh.p(v[3]);
        wrong += orient3d(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y,
                          d.z) != 1 ||
                 mesh.volume(t) <= 0.0;
    }
    check(wrong == 0,
          "orient3d is positive on every tet TetMesh orients positively");
}

} // namespace

int main()
{
    test_random();
    test_near_plane();
    test_on_plane();
    test_symmetry();
    test_tetmesh_convention();
    return report_checks();
}
