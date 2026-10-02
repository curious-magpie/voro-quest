// side1 (kernel/side.h): is a domain vertex q on seed i's side of the bisector
// of i and m? The sign of f_im(q) = pi_i(q) - pi_m(q), the difference of power
// distances pi_a(x) = |x - p_a|^2 - w_a: -1 when q is on i's side (closer to i
// in power distance), +1 on m's side, 0 exactly on the bisector. Degree 2.
//
// Against the oracle of side_support.h (the power distances from their
// definition, in rationals), on the driver's input families -- random, a
// coarse lattice full of ties, the lattice far from the origin, ties nudged by
// an ulp -- plus a few answers worked by hand, and antisymmetry: swapping i
// and m flips the sign.
#include "side_support.h"

#include "kernel/side.h"

namespace
{

int call(const Case &c)
{
    return side1(c.s[0].p, c.s[0].w, c.s[1].p, c.s[1].w, c.q[0]);
}

void test_by_hand()
{
    const glm::dvec3 o(0, 0, 0), m(2, 0, 0);
    check(side1(o, 0, m, 0, {0.5, 0, 0}) == -1, "side1: near i, -1");
    check(side1(o, 0, m, 0, {1.5, 0, 0}) == 1, "side1: near m, +1");
    check(side1(o, 0, m, 0, {1, 5, -7}) == 0, "side1: on the bisector, 0");
    // A weight pulls the bisector towards the other seed: with w_m = 4,
    // pi_m(1, 0, 0) = 1 - 4 < 1 = pi_i(1, 0, 0).
    check(side1(o, 0, m, 4, {1, 0, 0}) == 1,
          "side1: a weight on m moves the midpoint to m's side");
}

void test_antisymmetry()
{
    SideTally t;
    std::mt19937_64 rng(101);
    for (int it = 0; it < 50000; ++it)
    {
        const Case c = make_case(rng, 2, 1, it % 2 == 0);
        const Case swapped{{c.s[1], c.s[0]}, c.q};
        t.add(call(swapped) == -call(c), c);
    }
    t.report("side1: swapping i and m flips the sign");
}

} // namespace

int main()
{
    test_by_hand();
    check_predicate("side1", 2, 1, call, side1_counts(), {0.99, 0.8}, 100);
    test_antisymmetry();
    return report_checks();
}
