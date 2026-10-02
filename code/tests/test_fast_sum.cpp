// fast_sum (kernel/expansion.h): Shewchuk's Fast-Expansion-Sum, against exact
// rationals and against sum(), the slow version step 2 established.
//
// fast_sum has to give the same *value* as sum -- not necessarily the same
// components, since the two algorithms carry differently -- and a well formed
// expansion. That much is Shewchuk's Lemma 16, for strongly nonoverlapping
// inputs, and every input here is one: grow(), sum() and scale() give
// nonadjacent results.
//
// What is *not* checked any more is Theorem 13's further promise, that the
// result is strongly nonoverlapping again: it is false
// (test_dynamic_expansion.cpp pins a counterexample). So a fast_sum result is
// not known to be a valid input to the next fast_sum, and the operators sum
// with linear_sum instead (test_linear_sum.cpp). fast_sum stays, as the
// algorithm the counterexample is about. The tree below feeds fast_sum its own
// results anyway; it is evidence, not a theorem.
//
// Inputs come from grow(), sum() and scale().
#include <cmath>
#include <random>
#include <string>

#include "expansion_support.h"

namespace
{

// Every property of one fast_sum result, against the exact value it should
// have and the value sum() gives.
struct Checks
{
    Tally exact, form, agrees;

    template <std::size_t M, std::size_t N>
    void run(const Expansion<M> &a, const Expansion<N> &b)
    {
        const Expansion<M + N> f = fast_sum(a, b);
        const std::string ex = describe(a) + " + " + describe(b);
        const mpq_class want = value(a) + value(b);
        exact.add(value(f) == want, ex);
        form.add(well_formed(f), ex);
        agrees.add(value(f) == value(sum(a, b)), ex);
    }

    void report(const std::string &what) const
    {
        exact.report(("fast_sum, " + what + ": exactly a + b").c_str());
        form.report(("fast_sum, " + what + ": well formed").c_str());
        agrees.report(("fast_sum, " + what + ": the value sum() gives").c_str());
    }
};

void test_grown()
{
    Checks c;
    std::mt19937_64 rng(20);
    mpq_class ignored;
    for (int i = 0; i < 20000; ++i)
    {
        c.run(build<7>(rng, ignored), build<5>(rng, ignored));
        c.run(build<1>(rng, ignored), build<9>(rng, ignored));
        c.run(build<9>(rng, ignored), build<1>(rng, ignored));
    }
    c.report("grown inputs");
}

void test_made_other_ways()
{
    // Products and sums have different shapes of components from grown
    // expansions -- scale's come in pairs -- so they are inputs too.
    Checks c;
    std::mt19937_64 rng(21);
    mpq_class ignored;
    for (int i = 0; i < 20000; ++i)
    {
        const Expansion<5> a = build<5>(rng, ignored);
        const Expansion<4> b = build<4>(rng, ignored);
        const Expansion<10> p = scale(a, random_double(rng, -100, 100));
        const Expansion<9> s = sum(a, b);
        c.run(p, s);
        c.run(p, negate(p)); // cancels to zero, exactly
        c.run(s, a);
    }
    c.report("scaled and summed inputs");
}

void test_edges()
{
    Checks c;
    std::mt19937_64 rng(22);
    mpq_class ignored;
    const Expansion<6> a = build<6>(rng, ignored);
    const Expansion<3> zero;

    c.run(zero, zero);
    c.run(zero, a);
    c.run(a, zero);
    c.run(a, a); // equal magnitudes all the way: the merge's ties
    c.run(Expansion<1>(1.0), Expansion<1>(1.0));
    c.run(Expansion<1>(1.0), Expansion<1>(-1.0));
    c.run(Expansion<1>(std::ldexp(1.0, 100)), Expansion<1>(1.0));
    c.report("edge cases");

    check(fast_sum(zero, zero).size() == 0, "fast_sum of zeros has no components");
    check(fast_sum(a, negate(a)).size() == 0,
          "fast_sum of a and -a has no components");
}

// Sixteen expansions summed as a balanced tree of fast_sums: every level's
// inputs are the previous level's outputs. Only the first level is covered by
// a theorem (Lemma 16); above it the inputs are fast_sum results, which
// Theorem 13 was meant to cover. Evidence, then: no failure here has been
// seen, but none is ruled out.
void test_tree()
{
    Tally exact, form;
    std::mt19937_64 rng(23);
    for (int i = 0; i < 2000; ++i)
    {
        mpq_class v[16];
        Expansion<4> leaf[16];
        mpq_class want = 0;
        for (int k = 0; k < 16; ++k)
        {
            leaf[k] = build<4>(rng, v[k]);
            want += v[k];
        }
        Expansion<8> l1[8];
        for (int k = 0; k < 8; ++k)
            l1[k] = fast_sum(leaf[2 * k], leaf[2 * k + 1]);
        Expansion<16> l2[4];
        for (int k = 0; k < 4; ++k)
            l2[k] = fast_sum(l1[2 * k], l1[2 * k + 1]);
        const Expansion<32> l3a = fast_sum(l2[0], l2[1]);
        const Expansion<32> l3b = fast_sum(l2[2], l2[3]);
        const Expansion<64> top = fast_sum(l3a, l3b);

        exact.add(value(top) == want, describe(top));
        form.add(well_formed(top), describe(top));
    }
    exact.report("evidence: fast_sum tree, exactly the sum of sixteen");
    form.report("evidence: fast_sum tree, well formed");
}

} // namespace

int main()
{
    test_grown();
    test_made_other_ways();
    test_edges();
    test_tree();
    return report_checks();
}
