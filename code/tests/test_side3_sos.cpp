// side3_sos (kernel/side.h): side3 under Simulation of Simplicity. Never 0;
// equal to side3 where side3 is not 0. The denominator does not change under
// the perturbation, and the numerator becomes linear in the eps (a row plus a
// multiple of the ones row; terms with two such rows vanish):
//
//   num + a det(1; J; K) + b det(M; 1; K) + c det(M; J; 1)
//
// with a = eps_m - eps_i, b = eps_j - eps_i, c = eps_k - eps_i, and M, J, K
// the rows of f_im, f_ij, f_ik over the corners. det(1; J; K) is the
// denominator, so eps_m's coefficient is never zero.
//
// Against the perturbed oracle on the shared driver, and the vertex's
// symmetries at ties: any order of the corners, any of i, j, k as i.
#include <algorithm>

#include "sos_support.h"

#include "kernel/side.h"

namespace
{

int unperturbed(const Case &c)
{
    return side3(c.s[0].p, c.s[0].w, c.s[1].p, c.s[1].w, c.s[2].p, c.s[2].w,
                 c.s[3].p, c.s[3].w, c.q[0], c.q[1], c.q[2]);
}

int sos(const Case &c)
{
    return side3_sos(c.s[0].id, c.s[0].p, c.s[0].w, c.s[1].id, c.s[1].p,
                     c.s[1].w, c.s[2].id, c.s[2].p, c.s[2].w, c.s[3].id,
                     c.s[3].p, c.s[3].w, c.q[0], c.q[1], c.q[2]);
}

void test_symmetry()
{
    SideTally corners, seeds;
    std::mt19937_64 rng(311);
    for (int it = 0; it < 20000; ++it)
    {
        const Case c = with_ids(make_case(rng, 4, 3, true), rng);
        if (!oracle(c))
            continue;
        const int s = sos(c);
        bool same_q = true, same_s = true;
        std::vector<int> order{0, 1, 2};
        do
            same_q &=
                sos(Case{c.s, {c.q[order[0]], c.q[order[1]], c.q[order[2]]}}) ==
                s;
        while (std::next_permutation(order.begin(), order.end()));
        order = {0, 1, 2};
        do
            same_s &=
                sos(Case{{c.s[order[0]], c.s[order[1]], c.s[order[2]], c.s[3]},
                         c.q}) == s;
        while (std::next_permutation(order.begin(), order.end()));
        corners.add(same_q, c);
        seeds.add(same_s, c);
    }
    corners.report("side3_sos: any order of the corners, ties included");
    seeds.report("side3_sos: any of i, j, k as i, ties included");
}

} // namespace

int main()
{
    check_perturbed("side3_sos", 4, 3, unperturbed, sos, side3_counts(), 310);
    test_symmetry();
    return report_checks();
}
