// side1_sos (kernel/side.h): side1 under Simulation of Simplicity. The same
// question -- is q on i's side of the bisector of i and m? -- asked of the
// imaginary input where every weight w_a is w_a + eps_a, eps_0 >> eps_1 >>
// ... by seed index. It never answers 0. Where side1 is not 0 it gives the
// same answer; at a tie, f_im(q) + eps_m - eps_i decides, so the seed with
// the smaller index wins: its extra weight is infinitely larger.
//
// Against the perturbed oracle of sos_support.h on the shared driver, by hand,
// and antisymmetry at ties: swapping i and m (with their indices) flips the
// answer, ties included.
#include "sos_support.h"

#include "kernel/side.h"

namespace
{

int unperturbed(const Case &c)
{
    return side1(c.s[0].p, c.s[0].w, c.s[1].p, c.s[1].w, c.q[0]);
}

int sos(const Case &c)
{
    return side1_sos(c.s[0].id, c.s[0].p, c.s[0].w, c.s[1].id, c.s[1].p,
                     c.s[1].w, c.q[0]);
}

void test_by_hand()
{
    const glm::dvec3 o(0, 0, 0), m(2, 0, 0), on(1, 5, -7);
    check(side1_sos(0, o, 0, 1, m, 0, on) == -1,
          "side1_sos: a tie goes to the smaller index (i = 0): -1");
    check(side1_sos(1, o, 0, 0, m, 0, on) == 1,
          "side1_sos: a tie goes to the smaller index (m = 0): +1");
    check(side1_sos(5, o, 0, 9, m, 0, {0.5, 0, 0}) == -1,
          "side1_sos: no tie, the exact answer");
}

void test_antisymmetry()
{
    SideTally t;
    std::mt19937_64 rng(111);
    for (int it = 0; it < 30000; ++it)
    {
        const Case c = with_ids(make_case(rng, 2, 1, true), rng);
        t.add(sos(Case{{c.s[1], c.s[0]}, c.q}) == -sos(c), c);
    }
    t.report("side1_sos: swapping i and m flips the sign, ties included");
}

} // namespace

int main()
{
    test_by_hand();
    check_perturbed("side1_sos", 2, 1, unperturbed, sos, side1_counts(), 110);
    test_antisymmetry();
    return report_checks();
}
