// Where seeds come from.
//
// Only the uniform sampler exists yet. The lattices, the samplers on a surface
// and the readers the design lists (section 4.1) arrive with the milestones
// that test against them.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <glm/glm.hpp>

#include "domain/box.h"

// `count` points uniformly distributed in `box`, drawn from `seed`.
//
// The same seed gives the same points, bit for bit, on every machine and with
// every standard library: the generator is std::mt19937_64, whose output the
// standard fixes, and the conversion to [0, 1) is done here rather than by
// std::uniform_real_distribution, whose algorithm each library chooses for
// itself. The build keeps the final multiply-add unfused (-ffp-contract=off
// in CMakeLists.txt), so it rounds the same way everywhere too. A diagram that
// goes wrong on one machine can then be rebuilt on another from two numbers.
std::vector<glm::dvec3> sample_uniform(const Box &box, size_t count,
                                       uint64_t seed);
