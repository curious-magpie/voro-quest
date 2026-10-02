// side4 (kernel/side.h): is the Voronoi vertex of seeds i, j, k, l -- the
// point with equal power distance to all four -- on i's side of the bisector
// of i and m? -1 on i's side, +1 on m's, 0 when the five seeds are on one
// power sphere. With u_a = 2 (p_a - p_i) and c_a = |p_a|^2 - |p_i|^2 + w_i -
// w_a (or the same in coordinates relative to p_i),
//
//   f_im(x) = - det4[(u_a, c_a) for a = j, k, l, m] / det3[u_j; u_k; u_l]
//                                                      (design, 6.2)
//
// degree 5 over degree 3. Precondition: det3 is not zero (i, j, k, l not
// coplanar, so the vertex exists).
//
// Against the oracle (x solved exactly from the three bisector equations,
// then the power distances) on the driver's inputs; answers by hand; and the
// vertex's symmetry: any order of i, j, k, l is the same point with the same
// power distance, so the same sign.
#include <algorithm>

#include "side_support.h"

#include "kernel/side.h"

namespace
{

int call(const Case &c)
{
    return side4(c.s[0].p, c.s[0].w, c.s[1].p, c.s[1].w, c.s[2].p, c.s[2].w,
                 c.s[3].p, c.s[3].w, c.s[4].p, c.s[4].w);
}

void test_by_hand()
{
    // The Voronoi vertex of the origin and the three unit points scaled by 2
    // is x = (1,1,1), where pi_i(x) = 3.
    const glm::dvec3 pi(0, 0, 0), pj(2, 0, 0), pk(0, 2, 0), pl(0, 0, 2);
    check(side4(pi, 0, pj, 0, pk, 0, pl, 0, {3, 3, 3}, 0) == -1,
          "side4: m farther from x than i, -1");
    check(side4(pi, 0, pj, 0, pk, 0, pl, 0, {1, 1, 2}, 0) == 1,
          "side4: m closer to x than i, +1");
    check(side4(pi, 0, pj, 0, pk, 0, pl, 0, {2, 2, 2}, 0) == 0,
          "side4: five seeds on one sphere, 0");
    // Weights: w_m = 3 makes pi_m(x) = 12 - 3 = 9 > 3, still -1; w_m = 10
    // makes it 2 < 3.
    check(side4(pi, 0, pj, 0, pk, 0, pl, 0, {3, 3, 3}, 10) == 1,
          "side4: a large weight on m takes the vertex to m's side");
}

void test_symmetry()
{
    SideTally t;
    std::mt19937_64 rng(401);
    for (int it = 0; it < 20000; ++it)
    {
        const Case c = make_case(rng, 5, 0, it % 2 == 0);
        if (!oracle(c))
            continue;
        const int s = call(c);
        std::vector<int> order{0, 1, 2, 3};
        bool same = true;
        while (std::next_permutation(order.begin(), order.end()))
            same &= call(Case{{c.s[order[0]], c.s[order[1]], c.s[order[2]],
                               c.s[order[3]], c.s[4]},
                              {}}) == s;
        t.add(same, c);
    }
    t.report("side4: any order of i, j, k, l, the same sign");
}

} // namespace

int main()
{
    test_by_hand();
    check_predicate("side4", 5, 0, call, side4_counts(), {0.99, 0.8}, 400);
    test_symmetry();
    return report_checks();
}
