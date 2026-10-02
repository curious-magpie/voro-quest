// side2 (kernel/side.h): where the domain edge q0 q1 crosses the bisector of
// seeds i and j, is that point on i's side of the bisector of i and m? -1 on
// i's side, +1 on m's, 0 on it. The predicate takes the sign of
//
//   f_im(x) = det[[f_im(q0), f_im(q1)], [f_ij(q0), f_ij(q1)]]
//           / det[[1, 1], [f_ij(q0), f_ij(q1)]]        (design, 6.2)
//
// degree 4 over degree 2, without ever constructing x. Precondition: the
// denominator is not zero, i.e. the edge is not parallel to the bisector of
// i and j (the vertex exists). The tests never ask about such an edge.
//
// Against the oracle (x constructed exactly, then the power distances), on
// the driver's inputs; answers by hand; and two symmetries: the edge read the
// other way round, and i and j exchanged -- at x, pi_i = pi_j, so f_im(x) and
// f_jm(x) are the same number.
#include "side_support.h"

#include "kernel/side.h"

namespace
{

int call(const Case &c)
{
    return side2(c.s[0].p, c.s[0].w, c.s[1].p, c.s[1].w, c.s[2].p, c.s[2].w,
                 c.q[0], c.q[1]);
}

void test_by_hand()
{
    // The edge (-1,0,0)-(3,0,0) crosses the bisector of (0,0,0) and (2,0,0)
    // at x = (1,0,0), where pi_i(x) = 1.
    const glm::dvec3 pi(0, 0, 0), pj(2, 0, 0), q0(-1, 0, 0), q1(3, 0, 0);
    check(side2(pi, 0, pj, 0, {5, 0, 0}, 0, q0, q1) == -1,
          "side2: m far away, -1");
    check(side2(pi, 0, pj, 0, {1, 0.5, 0}, 0, q0, q1) == 1,
          "side2: m right next to x, +1");
    check(side2(pi, 0, pj, 0, {1, 1, 0}, 0, q0, q1) == 0,
          "side2: m as far from x as i is, 0");
    check(side2(pi, 0, pj, 0, {5, 0, 0}, 0, q1, q0) == -1,
          "side2: the edge read backwards, the same answer");
}

void test_symmetry()
{
    SideTally edge, ij;
    std::mt19937_64 rng(201);
    for (int it = 0; it < 50000; ++it)
    {
        const Case c = make_case(rng, 3, 2, it % 2 == 0);
        if (!oracle(c))
            continue;
        const int s = call(c);
        edge.add(call(Case{c.s, {c.q[1], c.q[0]}}) == s, c);
        ij.add(call(Case{{c.s[1], c.s[0], c.s[2]}, c.q}) == s, c);
    }
    edge.report("side2: the edge read the other way, the same sign");
    ij.report("side2: i and j exchanged, the same sign");
}

} // namespace

int main()
{
    test_by_hand();
    check_predicate("side2", 3, 2, call, side2_counts(), {0.99, 0.8}, 200);
    test_symmetry();
    return report_checks();
}
