// side2_sos (kernel/side.h): side2 under Simulation of Simplicity. Never 0;
// equal to side2 where side2 is not 0. At a tie the numerator becomes
//
//   num + eps_m (j1 - j0) + eps_j (m0 - m1) - eps_i ((j1 - j0) + (m0 - m1))
//
// (the denominator j1 - j0 does not change: its eps cancel), and the first
// nonzero coefficient in increasing seed index decides; eps_m's coefficient
// is the denominator, so there always is one.
//
// Against the perturbed oracle on the shared driver, and the vertex's
// symmetries at ties: the edge read backwards, and i and j exchanged (on the
// perturbed input too, x is where pi_i = pi_j).
#include "sos_support.h"

#include "kernel/side.h"

namespace
{

int unperturbed(const Case &c)
{
    return side2(c.s[0].p, c.s[0].w, c.s[1].p, c.s[1].w, c.s[2].p, c.s[2].w,
                 c.q[0], c.q[1]);
}

int sos(const Case &c)
{
    return side2_sos(c.s[0].id, c.s[0].p, c.s[0].w, c.s[1].id, c.s[1].p,
                     c.s[1].w, c.s[2].id, c.s[2].p, c.s[2].w, c.q[0], c.q[1]);
}

void test_symmetry()
{
    SideTally edge, ij;
    std::mt19937_64 rng(211);
    for (int it = 0; it < 30000; ++it)
    {
        const Case c = with_ids(make_case(rng, 3, 2, true), rng);
        if (!oracle(c))
            continue;
        const int s = sos(c);
        edge.add(sos(Case{c.s, {c.q[1], c.q[0]}}) == s, c);
        ij.add(sos(Case{{c.s[1], c.s[0], c.s[2]}, c.q}) == s, c);
    }
    edge.report("side2_sos: the edge read the other way, ties included");
    ij.report("side2_sos: i and j exchanged, ties included");
}

} // namespace

int main()
{
    check_perturbed("side2_sos", 3, 2, unperturbed, sos, side2_counts(), 210);
    test_symmetry();
    return report_checks();
}
