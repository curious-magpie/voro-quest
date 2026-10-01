// The seeds: the point set's revision rules, the uniform sampler, and the
// session that ties them to the box.
//
// The sampler is held to its header's promise -- the same seed gives the same
// points, bit for bit, anywhere -- because a diagram that goes wrong has to be
// rebuildable from two numbers.
#include <cmath>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "test_support.h"

#include "points/point_set.h"
#include "points/samplers.h"
#include "session/rvd_session.h"

namespace
{

bool inside(const Box &box, const glm::dvec3 &p)
{
    for (int a = 0; a < 3; ++a)
        if (!(p[a] >= box.lo()[a] && p[a] <= box.hi()[a]))
            return false;
    return true;
}

bool same_bits(const std::vector<glm::dvec3> &a,
               const std::vector<glm::dvec3> &b)
{
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(),
                                     a.size() * sizeof(glm::dvec3)) == 0);
}

void test_point_set()
{
    PointSet s;
    check(s.empty() && s.revision() == 0, "an empty set has revision 0");

    s.assign({{0, 0, 0}, {1, 2, 3}, {4, 5, 6}});
    check(s.size() == 3 && s.revision() == 1,
          "assign replaces the set and bumps the revision");
    check(s.weights().size() == 3 && s.weights()[0] == 0.0 &&
              s.weights()[2] == 0.0,
          "assign gives every seed a zero weight");

    s.positions_mut()[1].x = 7.0;
    check(s.revision() == 2 && s.p(1).x == 7.0,
          "writing positions bumps the revision");
    s.weights_mut()[0] = 0.5;
    check(s.revision() == 3, "writing weights bumps the revision");
}

void test_sampler()
{
    Box box;
    std::string err;
    check(box.build({-2.0, 0.0, 10.0}, {1.0, 0.5, 14.0}, err),
          "the sample box builds");

    const size_t n = 100000;
    const std::vector<glm::dvec3> a = sample_uniform(box, n, 42);
    check(a.size() == n, "the sampler draws exactly the count asked for");

    bool all_inside = true;
    for (const glm::dvec3 &p : a)
        all_inside = all_inside && inside(box, p);
    check(all_inside, "every sample lies in the box");

    check(same_bits(a, sample_uniform(box, n, 42)),
          "the same seed gives the same bits");
    check(!same_bits(a, sample_uniform(box, n, 43)),
          "a different seed gives other points");

    const std::vector<glm::dvec3> head = sample_uniform(box, 100, 42);
    check(same_bits(head, std::vector<glm::dvec3>(a.begin(), a.begin() + 100)),
          "a smaller draw is the start of a larger one");
    check(sample_uniform(box, 0, 42).empty(), "a count of zero draws nothing");

    // Moments. The standard error of the mean is size / sqrt(12 n), about 1e-3
    // of the size here, so a 1e-2 tolerance is ten of them: a failure is a bug,
    // not bad luck. The variance of a uniform variable is size^2 / 12.
    glm::dvec3 mean(0.0), var(0.0);
    for (const glm::dvec3 &p : a)
        mean += p;
    mean /= double(n);
    for (const glm::dvec3 &p : a)
        var += (p - mean) * (p - mean);
    var /= double(n - 1);
    const glm::dvec3 size = box.hi() - box.lo();
    bool mean_ok = true, var_ok = true;
    for (int k = 0; k < 3; ++k)
    {
        mean_ok =
            mean_ok && std::fabs(mean[k] - box.center()[k]) < 1e-2 * size[k];
        const double want = size[k] * size[k] / 12.0;
        var_ok = var_ok && std::fabs(var[k] - want) < 2e-2 * want;
    }
    check(mean_ok, "the samples' mean is the box's centre");
    check(var_ok, "the samples' variance is a uniform one");

    // What "the same on every machine" rests on: the standard fixes this
    // generator's output ([rand.predef]: the 10000th draw of a default-seeded
    // mt19937_64). If this fails, the library is not the standard one.
    std::mt19937_64 rng;
    rng.discard(9999);
    check(rng() == 9981545732273789042ull,
          "mt19937_64 is the standard's generator");
}

void test_session()
{
    RvdSession session;
    check(session.box().volume() == 1.0, "a new session works in the unit box");
    check(session.seeds().size() == size_t(session.seed_params.count),
          "a new session has its default seeds");

    const uint64_t before = session.seeds().revision();
    session.seed_params.count = 10;
    session.regenerate_seeds();
    check(session.seeds().size() == 10 && session.seeds().revision() > before,
          "regenerating draws the count asked for, and bumps the revision");

    session.seed_params.count = -5;
    session.regenerate_seeds();
    check(session.seeds().empty(), "a negative count draws nothing");
    session.seed_params.count = 10;
    session.regenerate_seeds();

    std::string err;
    const uint64_t kept = session.seeds().revision();
    check(!session.set_box({0, 0, 0}, {0, 1, 1}, err) && !err.empty(),
          "the session refuses a flat box");
    check(session.seeds().revision() == kept && session.box().volume() == 1.0,
          "a refused box changes neither the box nor the seeds");

    check(session.set_box({5, 5, 5}, {6, 7, 8}, err),
          "the session takes a new box", err.c_str());
    bool all_inside = session.seeds().size() == 10;
    for (const glm::dvec3 &p : session.seeds().positions())
        all_inside = all_inside && inside(session.box(), p);
    check(all_inside, "a new box redraws the seeds inside it");
}

} // namespace

int main()
{
    test_point_set();
    test_sampler();
    test_session();
    return report_checks();
}
