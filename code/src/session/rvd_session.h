// Everything RVD is working on: the domain, the seeds, and the parameters that
// made them. The counterpart of mesh-dev's Document (design, section 8.1).
//
// Core side, with no GL and no viewer: a test, a script or the app drives it
// the same way. The app shows it through the adapters in app/, which read it
// once a frame and never write to it.
//
// In milestone M0 the domain is always a box and the seeds are always uniform.
// The current generation -- the frozen result of a diagram (design, section
// 5.9) -- joins them here with the first diagram, in M2.
#pragma once

#include <cstdint>
#include <string>

#include "domain/box.h"
#include "points/point_set.h"

class RvdSession
{
  public:
    // How the seeds are drawn. Edited freely by the Diagram window: nothing
    // happens until regenerate_seeds() is called.
    struct SeedParams
    {
        int count = 2000;
        uint64_t random_seed = 1;
    };

    // The unit box, and seeds drawn with the default parameters.
    RvdSession();

    const Box &box() const
    {
        return box_;
    }
    // Replaces the box and redraws the seeds inside it. False with a reason in
    // `err`, and nothing changed, if the box is not a box.
    bool set_box(const glm::dvec3 &lo, const glm::dvec3 &hi, std::string &err);

    const PointSet &seeds() const
    {
        return seeds_;
    }

    SeedParams seed_params;

    // Draws seed_params.count uniform seeds in the box from
    // seed_params.random_seed. A negative count draws none.
    void regenerate_seeds();

  private:
    Box box_;
    PointSet seeds_;
};
