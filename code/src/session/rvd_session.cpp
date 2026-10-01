#include "session/rvd_session.h"

#include <algorithm>

#include "points/samplers.h"

RvdSession::RvdSession()
{
    // The unit box cannot fail to build, so the error is not looked at.
    std::string err;
    box_.build(glm::dvec3(0.0), glm::dvec3(1.0), err);
    regenerate_seeds();
}

bool RvdSession::set_box(const glm::dvec3 &lo, const glm::dvec3 &hi,
                         std::string &err)
{
    // Built aside first, so a rejected box leaves the old one, and its seeds,
    // exactly as they were.
    Box box = box_;
    if (!box.build(lo, hi, err))
        return false;
    box_ = box;
    regenerate_seeds();
    return true;
}

void RvdSession::regenerate_seeds()
{
    const size_t count = size_t(std::max(seed_params.count, 0));
    seeds_.assign(sample_uniform(box_, count, seed_params.random_seed));
}
