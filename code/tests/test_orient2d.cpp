// orient2d (kernel/orient2d.h): the first predicate, against exact rationals.
//
// orient2d(a, b, c) is the sign of the determinant
//
//   | ax - cx   ay - cy |
//   | bx - cx   by - cy |
//
// positive when a, b, c turn counterclockwise, negative clockwise, zero when
// they are collinear. The oracle here is written independently, as the cross
// product (b - a) x (c - a) in rationals, which is the same number by algebra
// but not by code.
//
// The hard inputs are the ones that matter: points a few ulps off a line,
// where the double formula gets the sign wrong -- the test checks that it does,
// so that passing here means something -- and points exactly on a line, where
// the only right answer is zero.
#include <cmath>
#include <random>
#include <string>

#include "test_support.h"

#include "exact.h"
#include "kernel/orient2d.h"

namespace
{

// The oracle: the sign of (b - a) x (c - a), exactly.
int oracle(double ax, double ay, double bx, double by, double cx, double cy)
{
    const mpq_class abx = exact(bx) - exact(ax), aby = exact(by) - exact(ay);
    const mpq_class acx = exact(cx) - exact(ax), acy = exact(cy) - exact(ay);
    return sgn(mpq_class(abx * acy - aby * acx));
}

struct Point
{
    double x, y;
};

std::string describe(const Point &a, const Point &b, const Point &c)
{
    return "(" + show(a.x) + ", " + show(a.y) + ") (" + show(b.x) + ", " +
           show(b.y) + ") (" + show(c.x) + ", " + show(c.y) + ")";
}

// Counts disagreements with the oracle, keeping the first.
struct Tally
{
    long failures = 0;
    std::string first;

    void add(bool ok, const Point &a, const Point &b, const Point &c)
    {
        if (!ok && failures++ == 0)
            first = describe(a, b, c);
    }
    void report(const char *what) const
    {
        check(failures == 0, what, failures ? first.c_str() : nullptr);
    }
};

int exact_sign(const Point &a, const Point &b, const Point &c)
{
    return orient2d(a.x, a.y, b.x, b.y, c.x, c.y);
}

int oracle_sign(const Point &a, const Point &b, const Point &c)
{
    return oracle(a.x, a.y, b.x, b.y, c.x, c.y);
}

// The same formula, evaluated in plain doubles: fast, and wrong near zero.
int double_sign(const Point &a, const Point &b, const Point &c)
{
    const double d = orient2d_det(a.x, a.y, b.x, b.y, c.x, c.y);
    return (d > 0) - (d < 0);
}

void test_random()
{
    Tally t;
    std::mt19937_64 rng(40);
    for (int i = 0; i < 200000; ++i)
    {
        Point p[3];
        for (Point &q : p)
            q = {random_double(rng, -30, 30), random_double(rng, -30, 30)};
        t.add(exact_sign(p[0], p[1], p[2]) == oracle_sign(p[0], p[1], p[2]),
              p[0], p[1], p[2]);
    }
    t.report("orient2d: random points, the exact sign");
}

// Kettner, Mehlhorn, Pion, Schirra and Yap, "Classroom examples of robustness
// problems in geometric computations" (2008): the 256 x 256 grid of points one
// ulp apart around (0.5, 0.5), against the line through (12, 12) and
// (24, 24). Most of the grid is within rounding of that line.
void test_near_line()
{
    Tally t;
    long double_wrong = 0;
    const Point b{12.0, 12.0}, c{24.0, 24.0};
    const double ulp = std::ldexp(1.0, -53); // the spacing just below 1
    for (int i = 0; i < 256; ++i)
        for (int j = 0; j < 256; ++j)
        {
            const Point a{0.5 + i * ulp, 0.5 + j * ulp};
            const int want = oracle_sign(a, b, c);
            t.add(exact_sign(a, b, c) == want, a, b, c);
            double_wrong += double_sign(a, b, c) != want;
        }
    t.report("orient2d: a grid one ulp apart near a line, the exact sign");
    check(double_wrong > 0,
          "the double formula gets some of that grid wrong (or the test is "
          "too easy)");
}

// Exactly on a line: a, b, and c = a + k (b - a) with integers, so that every
// coordinate is an exact double and the determinant is exactly zero.
//
// Unlike the grid above, this is not where doubles go wrong: when the
// coordinate differences are exact, the determinant's two products are the
// same real number, round to the same double, and cancel to an exact zero.
// What this checks is the exact path's own zero -- every component cancelling,
// to an empty expansion -- on products large enough (up to about 2^94) to need
// several components on the way.
void test_on_line()
{
    Tally t;
    std::mt19937_64 rng(41);
    std::uniform_int_distribution<long> coord(-(1L << 40), 1L << 40);
    std::uniform_int_distribution<int> k(-50, 50);
    for (int i = 0; i < 100000; ++i)
    {
        const Point a{double(coord(rng)), double(coord(rng))};
        const Point b{double(coord(rng)), double(coord(rng))};
        const int s = k(rng);
        const Point c{a.x + s * (b.x - a.x), a.y + s * (b.y - a.y)};
        // Scaled by a power of two to leave the integers, exactly.
        const double f = std::ldexp(1.0, int(rng() % 120) - 60);
        const Point fa{a.x * f, a.y * f}, fb{b.x * f, b.y * f}, fc{c.x * f, c.y * f};
        t.add(exact_sign(fa, fb, fc) == 0, fa, fb, fc);
    }
    t.report("orient2d: exactly collinear points give zero");
}

// Symmetries a sign must have: swapping two points flips it, rotating the
// three keeps it. Checked far from the origin, where cancellation is worst.
void test_symmetry()
{
    Tally swap, rotate, far;
    std::mt19937_64 rng(42);
    for (int i = 0; i < 100000; ++i)
    {
        Point p[3];
        for (Point &q : p)
            q = {1e6 + random_double(rng, -20, -10),
                 1e6 + random_double(rng, -20, -10)};
        const int s = exact_sign(p[0], p[1], p[2]);
        swap.add(exact_sign(p[1], p[0], p[2]) == -s, p[0], p[1], p[2]);
        rotate.add(exact_sign(p[1], p[2], p[0]) == s, p[0], p[1], p[2]);
        far.add(s == oracle_sign(p[0], p[1], p[2]), p[0], p[1], p[2]);
    }
    swap.report("orient2d: swapping two points flips the sign");
    rotate.report("orient2d: rotating the points keeps the sign");
    far.report("orient2d: points near 1e6, the exact sign");
}

} // namespace

int main()
{
    test_random();
    test_near_line();
    test_on_line();
    test_symmetry();
    return report_checks();
}
