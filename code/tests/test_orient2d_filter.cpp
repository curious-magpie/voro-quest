// The filtered orient2d (kernel/orient2d.h): the same exact answers, mostly
// without the exact path.
//
// test_orient2d.cpp already holds orient2d to the oracle on every hard input.
// This suite checks what the filter adds, through the predicate's counters:
//
//   on ordinary input -- random points, nowhere near collinear -- the filter
//   decides almost every call, which is the whole point of having one
//
//   on hard input -- a grid one ulp from a line, exactly collinear points --
//   calls do reach the exact path, and the answers are still the oracle's
//
//   every call is counted exactly once, as filtered or exact
#include <cmath>
#include <random>

#include "test_support.h"

#include "exact.h"
#include "kernel/orient2d.h"

namespace
{

int oracle(double ax, double ay, double bx, double by, double cx, double cy)
{
    const mpq_class abx = exact(bx) - exact(ax), aby = exact(by) - exact(ay);
    const mpq_class acx = exact(cx) - exact(ax), acy = exact(cy) - exact(ay);
    return sgn(mpq_class(abx * acy - aby * acx));
}

void test_random()
{
    orient2d_counts() = PredicateCounts{};
    std::mt19937_64 rng(60);
    const long n = 200000;
    long wrong = 0;
    for (long i = 0; i < n; ++i)
    {
        double p[6];
        for (double &v : p)
            v = random_double(rng, -30, 30);
        wrong += orient2d(p[0], p[1], p[2], p[3], p[4], p[5]) !=
                 oracle(p[0], p[1], p[2], p[3], p[4], p[5]);
    }
    const PredicateCounts c = orient2d_counts();
    check(wrong == 0, "filtered orient2d: random points, the exact sign");
    check(c.filtered + c.exact == uint64_t(n),
          "every call is counted once, as filtered or exact");
    check(c.filtered >= uint64_t(0.999 * n),
          "the filter decides at least 99.9% of random points");
}

void test_hard()
{
    // The grid of test_orient2d.cpp, one ulp apart near the line through
    // (12, 12) and (24, 24): where doubles fail, so the filter must not
    // decide everything.
    orient2d_counts() = PredicateCounts{};
    const double ulp = std::ldexp(1.0, -53);
    long wrong = 0;
    for (int i = 0; i < 256; ++i)
        for (int j = 0; j < 256; ++j)
        {
            const double ax = 0.5 + i * ulp, ay = 0.5 + j * ulp;
            wrong += orient2d(ax, ay, 12, 12, 24, 24) !=
                     oracle(ax, ay, 12, 12, 24, 24);
        }
    const PredicateCounts grid = orient2d_counts();
    check(wrong == 0, "filtered orient2d: the near-line grid, the exact sign");
    check(grid.exact > 0, "the near-line grid reaches the exact path");

    // Exactly collinear: the determinant is zero and its bound is not, so
    // no filter can call it, and every one must go exact.
    orient2d_counts() = PredicateCounts{};
    std::mt19937_64 rng(61);
    std::uniform_int_distribution<long> coord(-(1L << 40), 1L << 40);
    long nonzero = 0;
    const long n = 20000;
    for (long i = 0; i < n; ++i)
    {
        const double ax = double(coord(rng)), ay = double(coord(rng));
        const double bx = double(coord(rng)), by = double(coord(rng));
        const double cx = ax + 3 * (bx - ax), cy = ay + 3 * (by - ay);
        nonzero += orient2d(ax, ay, bx, by, cx, cy) != 0;
    }
    const PredicateCounts line = orient2d_counts();
    check(nonzero == 0, "filtered orient2d: collinear points give zero");
    check(line.exact == uint64_t(n),
          "collinear points always reach the exact path");
}

} // namespace

int main()
{
    test_random();
    test_hard();
    return report_checks();
}
