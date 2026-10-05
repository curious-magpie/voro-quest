// side4_sos (kernel/side.h): side4 under Simulation of Simplicity. Never 0;
// equal to side4 where side4 is not 0. Each c_a = |d_a|^2 + w_i - w_a gains
// eps_i - eps_a; det4 is linear in its c column and det3 holds no c, so
//
//   num = -det4 + sum over a in {j, k, l, m} of eps_a C_a - eps_i (sum C_a)
//
// with C_a the cofactor of c_a in det4 (signs - + - + down the column). C_m
// is det3, the denominator: never zero, so the chain always ends.
//
// Against the perturbed oracle on the shared driver, and the vertex's
// symmetry at ties: any order of i, j, k, l (with their indices) is the same
// point of the same perturbed input.
#include <algorithm>

#include "sos_support.h"

#include "kernel/side.h"

namespace
{

int unperturbed(const Case &c)
{
    return side4(c.s[0].p, c.s[0].w, c.s[1].p, c.s[1].w, c.s[2].p, c.s[2].w,
                 c.s[3].p, c.s[3].w, c.s[4].p, c.s[4].w);
}

int sos(const Case &c)
{
    return side4_sos(c.s[0].id, c.s[0].p, c.s[0].w, c.s[1].id, c.s[1].p,
                     c.s[1].w, c.s[2].id, c.s[2].p, c.s[2].w, c.s[3].id,
                     c.s[3].p, c.s[3].w, c.s[4].id, c.s[4].p, c.s[4].w);
}

void test_symmetry()
{
    SideTally t;
    std::mt19937_64 rng(411);
    for (int it = 0; it < 10000; ++it)
    {
        const Case c = with_ids(make_case(rng, 5, 0, true), rng);
        if (!oracle(c))
            continue;
        const int s = sos(c);
        std::vector<int> order{0, 1, 2, 3};
        bool same = true;
        while (std::next_permutation(order.begin(), order.end()))
            same &= sos(Case{{c.s[order[0]], c.s[order[1]], c.s[order[2]],
                              c.s[order[3]], c.s[4]},
                             {}}) == s;
        t.add(same, c);
    }
    t.report("side4_sos: any order of i, j, k, l, ties included");
}

} // namespace

int main()
{
    check_perturbed("side4_sos", 5, 0, unperturbed, sos, side4_counts(), 410);
    test_symmetry();
    return report_checks();
}
