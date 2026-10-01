// The seeds a diagram is built from: a position and a weight each.
//
// A seed *is* its index. There is no id and no attribute store: whatever
// belongs to a seed -- its cell's mass, its colour class, its neighbours -- is
// a plain array indexed the same way, owned by whoever computes it.
//
// Three decisions, each cheaper now than later:
//
//   Double, never rounded, never recentred. These are the input numbers every
//   predicate is evaluated from (design, section 6.1), so they are kept exactly
//   as given. The float32 recentred copy is the viewer's, and only the
//   viewer's.
//
//   A weight per seed from the start, zero for a Voronoi diagram. Power
//   diagrams and optimal transport need them, and adding them later would mean
//   touching every predicate (design, section 3.5).
//
//   A revision, bumped on every write and never zero once anything has been
//   assigned, so a cache that has never seen this set cannot appear up to date.
//   The viewer's layer uses it, and so will every structure built from the
//   seeds.
#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include <glm/glm.hpp>

// The positions are handed out as one flat array: to the viewer as a span of
// dvec3, and later to numpy as an n x 3 array. Both assume three packed
// doubles, so say so out loud, as tet_mesh.h does.
static_assert(sizeof(glm::dvec3) == 3 * sizeof(double),
              "PointSet assumes dvec3 is three packed doubles");

class PointSet
{
  public:
    size_t size() const
    {
        return pos_.size();
    }
    bool empty() const
    {
        return pos_.empty();
    }

    const glm::dvec3 &p(size_t i) const
    {
        return pos_[i];
    }
    const std::vector<glm::dvec3> &positions() const
    {
        return pos_;
    }
    // One per seed, and all zero until something sets them: a plain Voronoi
    // diagram is a power diagram whose weights are equal.
    const std::vector<double> &weights() const
    {
        return weights_;
    }

    // Replaces the whole set, weights reset to zero.
    void assign(std::vector<glm::dvec3> positions)
    {
        pos_ = std::move(positions);
        weights_.assign(pos_.size(), 0.0);
        ++revision_;
    }

    // Seeds move (a Lloyd step does exactly this); their number does not.
    // Writing through either of these bumps the revision.
    std::vector<glm::dvec3> &positions_mut()
    {
        ++revision_;
        return pos_;
    }
    std::vector<double> &weights_mut()
    {
        ++revision_;
        return weights_;
    }

    uint64_t revision() const
    {
        return revision_;
    }

  private:
    std::vector<glm::dvec3> pos_;
    std::vector<double> weights_; // one per seed
    uint64_t revision_ = 0;
};
