#include "domain/box.h"

#include <cmath>

// Two triangles per face, each face listed as its outward axis: -x, +x, -y,
// +y, -z, +z. The winding is what makes the viewer's backface culling and its
// stencil cap work on the box, and tests/test_domain.cpp checks every normal.
const uint32_t Box::kTriangles[12][3] = {
    {0, 4, 6}, {0, 6, 2}, {1, 3, 7}, {1, 7, 5}, {0, 1, 5}, {0, 5, 4},
    {2, 6, 7}, {2, 7, 3}, {0, 2, 3}, {0, 3, 1}, {4, 5, 7}, {4, 7, 6}};

bool Box::build(const glm::dvec3 &lo, const glm::dvec3 &hi, std::string &err)
{
    for (int a = 0; a < 3; ++a)
    {
        if (!std::isfinite(lo[a]) || !std::isfinite(hi[a]))
        {
            err = "box corners must be finite";
            return false;
        }
        // Written as a negation so that equal corners -- a flat box -- fail
        // too.
        if (!(lo[a] < hi[a]))
        {
            err = "box is flat or inside out along axis " + std::to_string(a);
            return false;
        }
    }

    lo_ = lo;
    hi_ = hi;
    for (int c = 0; c < kCorners; ++c)
        corners_[size_t(c)] =
            glm::dvec3((c & 1) ? hi.x : lo.x, (c & 2) ? hi.y : lo.y,
                       (c & 4) ? hi.z : lo.z);
    for (size_t t = 0; t < 12; ++t)
        for (size_t k = 0; k < 3; ++k)
            triangles_[3 * t + k] = kTriangles[t][k];

    ++revision_;
    return true;
}

double Box::volume() const
{
    const glm::dvec3 d = hi_ - lo_;
    return d.x * d.y * d.z;
}
