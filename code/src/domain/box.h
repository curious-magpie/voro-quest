// The simplest domain: an axis-aligned box.
//
// It is the domain of a plain Voronoi diagram, which is the whole point of
// having it as a domain rather than as a special case (design, section 11.2,
// rule 1): a Voronoi diagram in a box is a restricted Voronoi diagram whose
// domain is one convex piece. Everything later milestones do to a triangle or
// a tet they will do to this box first.
//
// Its planes are axis-aligned and its corners are the input doubles
// themselves, so every coefficient a predicate reads off it is exact.
//
// What is here now is what milestone M0 needs: the corners, and the twelve
// triangles of its boundary to draw. The features the clippers name vertices by
// -- corners, edges and faces with ids -- come with the first clipper, in M2.
#pragma once

#include <array>
#include <cstdint>
#include <string>

#include <glm/glm.hpp>

class Box
{
  public:
    // Corner c has bit 1 set for the high x, bit 2 for the high y, bit 4 for
    // the high z -- the numbering mesh-dev's test cubes use.
    static constexpr int kCorners = 8;

    // The boundary as twelve triangles over the corners, each wound so that its
    // normal points out of the box.
    static const uint32_t kTriangles[12][3];

    // The box from `lo` to `hi`. False with a reason in `err` if the corners
    // are not finite, or if the box is flat or inside out on any axis: a domain
    // with no volume has no diagram to restrict to.
    bool build(const glm::dvec3 &lo, const glm::dvec3 &hi, std::string &err);

    const glm::dvec3 &lo() const
    {
        return lo_;
    }
    const glm::dvec3 &hi() const
    {
        return hi_;
    }
    glm::dvec3 center() const
    {
        return (lo_ + hi_) * 0.5;
    }
    double volume() const;

    // The eight corners, as positions to draw. Rebuilt with the box.
    const std::array<glm::dvec3, kCorners> &corners() const
    {
        return corners_;
    }
    // kTriangles flattened, three indices into corners() each, for the viewer.
    const std::array<uint32_t, 36> &boundary_triangles() const
    {
        return triangles_;
    }

    // Changes whenever build() succeeds, and never zero afterwards, so a cache
    // that has never seen this box cannot appear up to date.
    uint64_t revision() const
    {
        return revision_;
    }

  private:
    glm::dvec3 lo_{0.0}, hi_{0.0};
    std::array<glm::dvec3, kCorners> corners_{};
    std::array<uint32_t, 36> triangles_{};
    uint64_t revision_ = 0;
};
