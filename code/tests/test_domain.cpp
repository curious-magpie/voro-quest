// The box domain: what it accepts, and that its boundary is a closed surface
// wound outward -- which the viewer's culling and cap rely on today, and every
// clipper will rely on from M2.
//
// The closedness is checked by mesh-dev's own SurfaceMesh, built over the box's
// corners and triangles: the same test mesh-dev applies to every surface it
// loads, and the first use of the code RVD takes from it.
#include <cmath>
#include <string>
#include <vector>

#include "test_support.h"

#include "domain/box.h"
#include "mesh/surface_mesh.h"

namespace
{

void test_accepts_and_rejects()
{
    Box box;
    std::string err;
    check(box.revision() == 0, "an unbuilt box has revision 0");

    check(box.build(glm::dvec3(0.0), glm::dvec3(1.0), err),
          "the unit box builds");
    check(box.revision() == 1, "a build bumps the revision");
    check(box.volume() == 1.0, "the unit box has volume 1");

    // Corner c is high in x, y, z for bits 1, 2, 4.
    bool corners_ok = true;
    for (int c = 0; c < Box::kCorners; ++c)
    {
        const glm::dvec3 want(c & 1 ? 1.0 : 0.0, c & 2 ? 1.0 : 0.0,
                              c & 4 ? 1.0 : 0.0);
        corners_ok = corners_ok && box.corners()[size_t(c)] == want;
    }
    check(corners_ok, "corners follow the bit numbering");

    // Each of these must fail, and must leave the box it had untouched.
    const double inf = INFINITY, nan = NAN;
    const struct
    {
        glm::dvec3 lo, hi;
        const char *what;
    } bad[] = {
        {{0, 0, 0}, {0, 1, 1}, "a box flat in x is rejected"},
        {{0, 2, 0}, {1, 1, 1}, "a box inside out in y is rejected"},
        {{0, 0, nan}, {1, 1, 1}, "a NaN corner is rejected"},
        {{0, 0, 0}, {1, inf, 1}, "an infinite corner is rejected"},
    };
    for (const auto &b : bad)
    {
        err.clear();
        check(!box.build(b.lo, b.hi, err) && !err.empty(), b.what);
    }
    check(box.revision() == 1 && box.hi() == glm::dvec3(1.0),
          "a rejected build changes nothing");
}

void test_boundary()
{
    // Deliberately not a cube, and not at the origin, so that a wrong corner in
    // the table cannot hide behind a symmetry.
    const glm::dvec3 lo(-1.0, 0.0, 5.0), hi(2.0, 3.0, 5.5);
    Box box;
    std::string err;
    check(box.build(lo, hi, err), "an off-centre box builds", err.c_str());

    const glm::dvec3 d = hi - lo;
    const glm::dvec3 centre = box.center();
    const auto &p = box.corners();
    const auto &tris = box.boundary_triangles();

    double area = 0.0, volume = 0.0;
    bool outward = true, axis_aligned = true;
    for (size_t t = 0; t < 12; ++t)
    {
        const glm::dvec3 &a = p[tris[3 * t]], &b = p[tris[3 * t + 1]],
                         &c = p[tris[3 * t + 2]];
        const glm::dvec3 n = glm::cross(b - a, c - a); // twice the area vector
        area += 0.5 * glm::length(n);
        // Signed volume of the tet from the origin: summed over a closed
        // surface wound outward, it is the enclosed volume, with a positive
        // sign.
        volume += glm::dot(a, glm::cross(b, c)) / 6.0;

        outward = outward && glm::dot(n, (a + b + c) / 3.0 - centre) > 0.0;
        int nonzero = 0;
        for (int k = 0; k < 3; ++k)
            nonzero += n[k] != 0.0;
        axis_aligned = axis_aligned && nonzero == 1;
    }
    check(outward, "every boundary triangle faces out of the box");
    check(axis_aligned, "every boundary triangle lies in a face of the box");

    const double want_area = 2.0 * (d.x * d.y + d.y * d.z + d.z * d.x);
    check(std::fabs(area - want_area) <= 1e-12 * want_area,
          "the boundary's area is the box's");
    check(std::fabs(volume - box.volume()) <= 1e-12 * box.volume(),
          "the boundary encloses the box's volume, positively");

    // mesh-dev's verdict on the same triangles.
    SurfaceMesh surface;
    std::vector<glm::dvec3> pos(p.begin(), p.end());
    std::vector<uint32_t> ids(tris.begin(), tris.end());
    check(surface.build(std::move(pos), std::move(ids), err),
          "SurfaceMesh accepts the boundary", err.c_str());
    check(surface.closed(), "SurfaceMesh finds the boundary closed");
    check(surface.oriented(),
          "SurfaceMesh finds the boundary consistently wound");
    check(surface.degenerate_triangles() == 0,
          "no boundary triangle is degenerate");
}

} // namespace

int main()
{
    test_accepts_and_rejects();
    test_boundary();
    return report_checks();
}
