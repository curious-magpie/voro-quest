#include "points/samplers.h"

#include <random>

namespace
{

// The top 53 bits of a 64-bit draw, as a double in [0, 1): every value it can
// return is exactly representable, and 1 is not one of them.
double unit_interval(std::mt19937_64 &rng)
{
    return double(rng() >> 11) * 0x1.0p-53;
}

} // namespace

std::vector<glm::dvec3> sample_uniform(const Box &box, size_t count,
                                       uint64_t seed)
{
    std::mt19937_64 rng(seed);
    const glm::dvec3 lo = box.lo();
    const glm::dvec3 size = box.hi() - box.lo();

    // x, then y, then z for each point in turn, so that the first k points of a
    // larger draw are the k points of a smaller one with the same seed.
    std::vector<glm::dvec3> points(count);
    for (glm::dvec3 &p : points)
    {
        const double u = unit_interval(rng);
        const double v = unit_interval(rng);
        const double w = unit_interval(rng);
        p = lo + glm::dvec3(u, v, w) * size;
    }
    return points;
}
