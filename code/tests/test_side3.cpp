// side3 (kernel/side.h): where the domain triangle q0 q1 q2 crosses the
// bisectors of i and j and of i and k, is that point on i's side of the
// bisector of i and m? -1 on i's side, +1 on m's, 0 on it. The sign of
//
//   f_im(x) = det[[f_im(q0), f_im(q1), f_im(q2)],
//                 [f_ij(q0), f_ij(q1), f_ij(q2)],
//                 [f_ik(q0), f_ik(q1), f_ik(q2)]]
//           / det[[1, 1, 1],
//                 [f_ij(q0), f_ij(q1), f_ij(q2)],
//                 [f_ik(q0), f_ik(q1), f_ik(q2)]]      (design, 6.2)
//
// degree 6 over degree 4. Precondition: the denominator is not zero (the two
// bisectors cross the triangle's plane in one point).
//
// Against the oracle on the driver's inputs; answers by hand; and the
// symmetries of the vertex: it is the same point whatever the order of the
// corners, and whichever of i, j, k is called i (at x all three power
// distances are equal).
#include <algorithm>

#include "side_support.h"

#include "kernel/side.h"

namespace
{

int call(const Case &c)
{
    return side3(c.s[0].p, c.s[0].w, c.s[1].p, c.s[1].w, c.s[2].p, c.s[2].w,
                 c.s[3].p, c.s[3].w, c.q[0], c.q[1], c.q[2]);
}

void test_by_hand()
{
    // In the plane z = 0, the bisectors of (0,0,0) with (2,0,0) and with
    // (0,2,0) meet at x = (1,1,0), where pi_i(x) = 2.
    const glm::dvec3 pi(0, 0, 0), pj(2, 0, 0), pk(0, 2, 0);
    const glm::dvec3 q0(-4, -4, 0), q1(8, -4, 0), q2(-4, 8, 0);
    check(side3(pi, 0, pj, 0, pk, 0, {1, 1, 2}, 0, q0, q1, q2) == -1,
          "side3: m farther from x than i, -1");
    check(side3(pi, 0, pj, 0, pk, 0, {1, 1, 1}, 0, q0, q1, q2) == 1,
          "side3: m closer to x than i, +1");
    check(side3(pi, 0, pj, 0, pk, 0, {2, 2, 0}, 0, q0, q1, q2) == 0,
          "side3: m as far from x as i, 0");
}

void test_symmetry()
{
    SideTally corners, seeds;
    std::mt19937_64 rng(301);
    for (int it = 0; it < 30000; ++it)
    {
        const Case c = make_case(rng, 4, 3, it % 2 == 0);
        if (!oracle(c))
            continue;
        const int s = call(c);
        std::vector<glm::dvec3> q = c.q;
        std::vector<Site> ijk{c.s[0], c.s[1], c.s[2]};
        bool same_q = true, same_s = true;
        std::sort(
            q.begin(), q.end(), [](auto &a, auto &b)
            { return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z); });
        do
            same_q &= call(Case{c.s, q}) == s;
        while (std::next_permutation(
            q.begin(), q.end(), [](auto &a, auto &b)
            { return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z); }));
        // The three rotations of (i, j, k), and one swap.
        for (int r = 0; r < 3; ++r)
        {
            std::rotate(ijk.begin(), ijk.begin() + 1, ijk.end());
            same_s &= call(Case{{ijk[0], ijk[1], ijk[2], c.s[3]}, c.q}) == s;
        }
        same_s &= call(Case{{c.s[0], c.s[2], c.s[1], c.s[3]}, c.q}) == s;
        corners.add(same_q, c);
        seeds.add(same_s, c);
    }
    corners.report("side3: any order of the corners, the same sign");
    seeds.report("side3: any of i, j, k as i, the same sign");
}

} // namespace

int main()
{
    test_by_hand();
    check_predicate("side3", 4, 3, call, side3_counts(), {0.95, 0.8}, 300);
    test_symmetry();
    return report_checks();
}
