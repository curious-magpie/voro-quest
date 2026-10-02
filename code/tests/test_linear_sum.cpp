// linear_sum (kernel/expansion.h): Shewchuk's Linear-Expansion-Sum (Theorem
// 24, his Appendix A), the sum the operators use, against exact rationals.
//
// Why it replaced fast_sum in the operators: fast_sum needs strongly
// nonoverlapping inputs, and Theorem 13's promise that it gives one back is
// false (test_dynamic_expansion.cpp pins a counterexample). So a chain of
// fast_sums -- every product, every a + b + c -- rested on a step nobody has
// proved. Theorem 24 asks only what every Expansion already is,
// nonoverlapping and in increasing magnitude, and gives the same back, under
// any tie-breaking rule. So every check here is a theorem, chains included:
//
//   exact      value(linear_sum(a, b)) == value(a) + value(b)
//   form       the result is well formed: nonzero, nonoverlapping, increasing
//   any input  grown, scaled and summed expansions, and dense ones (touching
//              components of many bits: valid, but not strongly
//              nonoverlapping), on which fast_sum is not promised to work
//   chains     a tree of linear_sums, and linear_sum after scale after
//              linear_sum: each step's input is the last step's output
//   pinned     exact components on three inputs: Shewchuk's own example of
//              fast_sum producing overlap, and two counterexamples to Theorem
//              13. A correct sum by another algorithm (sum(), or
//              compress(fast_sum())) has the right value there but other
//              components, so these also check that this is Shewchuk's
//              algorithm
//   switch     a + b, a - b and product() are linear_sum now. These checks
//              fail until the operators and product() are switched over
#include <cmath>
#include <initializer_list>
#include <random>
#include <string>
#include <type_traits>

#include "expansion_support.h"

namespace
{

// An expansion from its components, smallest first.
template <std::size_t N>
Expansion<N> components(std::initializer_list<double> cs)
{
    Expansion<N> e;
    for (double c : cs)
        e.append(c);
    return e;
}

// Every property of one linear_sum result: its value, its form, its sign, and
// that sum() (Expansion-Sum, which also needs only nonoverlapping inputs)
// gets the same value.
struct Checks
{
    Tally exact, form, signs, agrees;

    template <std::size_t M, std::size_t N>
    void run(const Expansion<M> &a, const Expansion<N> &b)
    {
        const Expansion<M + N> s = linear_sum(a, b);
        const std::string ex = describe(a) + " + " + describe(b);
        const mpq_class want = value(a) + value(b);
        exact.add(value(s) == want, ex);
        form.add(well_formed(s), ex);
        signs.add(sign(s) == sgn(want), ex);
        agrees.add(value(s) == value(sum(a, b)), ex);
    }

    void report(const std::string &what) const
    {
        exact.report(("linear_sum, " + what + ": exactly a + b").c_str());
        form.report(("linear_sum, " + what + ": well formed").c_str());
        signs.report(("linear_sum, " + what + ": the right sign").c_str());
        agrees.report(
            ("linear_sum, " + what + ": the value sum() gives").c_str());
    }
};

void test_grown()
{
    Checks c;
    std::mt19937_64 rng(90);
    mpq_class ignored;
    for (int i = 0; i < 20000; ++i)
    {
        c.run(build<7>(rng, ignored), build<5>(rng, ignored));
        c.run(build<1>(rng, ignored), build<9>(rng, ignored));
        c.run(build<9>(rng, ignored), build<1>(rng, ignored));
        c.run(build<1>(rng, ignored), build<1>(rng, ignored));
        c.run(build<2>(rng, ignored), build<1>(rng, ignored));
    }
    c.report("grown inputs");
}

void test_made_other_ways()
{
    // Scaled results come in pairs of components, sums in other shapes.
    Checks c;
    std::mt19937_64 rng(91);
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
        c.run(fast_sum(a, b), p); // a fast_sum result is a valid input too
    }
    c.report("scaled and summed inputs");
}

void test_dense()
{
    // The inputs that set linear_sum apart: nonoverlapping, but touching.
    Checks c;
    std::mt19937_64 rng(92);
    mpq_class ignored;
    for (int i = 0; i < 40000; ++i)
    {
        const Expansion<6> a = dense<6>(rng, ignored);
        const Expansion<4> b = dense<4>(rng, ignored);
        c.run(a, b);
        c.run(b, a);
        c.run(a, a); // equal magnitudes all the way: the merge's ties
        c.run(a, negate(a));
        c.run(a, Expansion<1>(-approximate(a))); // cancels the top bits
        c.run(a, build<3>(rng, ignored));
        // Few-bit components, so that carries line up across both inputs.
        c.run(dense<8>(rng, ignored, 4), dense<8>(rng, ignored, 4));
    }
    c.report("dense inputs");
}

void test_edges()
{
    Checks c;
    std::mt19937_64 rng(93);
    mpq_class ignored;
    const Expansion<6> a = dense<6>(rng, ignored);
    const Expansion<3> zero;

    c.run(zero, zero);
    c.run(zero, a);
    c.run(a, zero);
    c.run(zero, Expansion<1>(3.0)); // one component in all: no pair to start
    c.run(Expansion<1>(3.0), zero);
    c.run(Expansion<1>(1.0), Expansion<1>(1.0));
    c.run(Expansion<1>(1.0), Expansion<1>(-1.0));
    c.run(Expansion<1>(std::ldexp(1.0, 100)), Expansion<1>(1.0));
    c.run(Expansion<1>(1.0), Expansion<1>(std::ldexp(1.0, -100)));
    c.report("edge cases");

    check(linear_sum(zero, zero).size() == 0,
          "linear_sum of zeros has no components");
    check(linear_sum(a, negate(a)).size() == 0,
          "linear_sum of a and -a has no components");
    check(same_components(linear_sum(zero, Expansion<1>(3.0)),
                          Expansion<1>(3.0)) &&
              same_components(linear_sum(Expansion<1>(3.0), zero),
                              Expansion<1>(3.0)),
          "linear_sum of one component and zero is that component");
    static_assert(
        std::is_same_v<decltype(linear_sum(Expansion<3>(), Expansion<2>())),
                       Expansion<5>>,
        "linear_sum has capacity M + N");
}

// Sixteen dense expansions summed as a balanced tree: every level's inputs
// are the previous level's outputs. With fast_sum this was the unproven
// case; with linear_sum each level is Theorem 24 again.
void test_tree()
{
    Tally exact_ok, form;
    std::mt19937_64 rng(94);
    for (int i = 0; i < 4000; ++i)
    {
        mpq_class v[16];
        Expansion<3> leaf[16];
        mpq_class want = 0;
        for (int k = 0; k < 16; ++k)
        {
            leaf[k] = (k % 2) ? dense<3>(rng, v[k]) : build<3>(rng, v[k]);
            want += v[k];
        }
        Expansion<6> l1[8];
        for (int k = 0; k < 8; ++k)
            l1[k] = linear_sum(leaf[2 * k], leaf[2 * k + 1]);
        Expansion<12> l2[4];
        for (int k = 0; k < 4; ++k)
            l2[k] = linear_sum(l1[2 * k], l1[2 * k + 1]);
        const Expansion<24> l3a = linear_sum(l2[0], l2[1]);
        const Expansion<24> l3b = linear_sum(l2[2], l2[3]);
        const Expansion<48> top = linear_sum(l3a, l3b);

        exact_ok.add(value(top) == want, describe(top));
        form.add(well_formed(top), describe(top));
    }
    exact_ok.report("linear_sum tree: exactly the sum of sixteen expansions");
    form.report("linear_sum tree: well formed");
}

// x -> x * d + c, three times, with linear_sum and scale alternating: what a
// predicate's formula does. Each step is a theorem (19 for scale, 24 here).
void test_chain()
{
    Tally exact_ok, form, signs;
    std::mt19937_64 rng(95);
    for (int i = 0; i < 20000; ++i)
    {
        mpq_class vx, vc;
        const Expansion<4> x0 = dense<4>(rng, vx, 20);
        const Expansion<3> c = dense<3>(rng, vc, 20);
        // Few-bit factors as well as full ones: they make carries line up.
        const double d = (i % 2) ? random_double(rng, -30, 30)
                                 : double(int(rng() % 31) - 15);
        const auto x1 = linear_sum(scale(x0, d), c);
        const auto x2 = linear_sum(scale(x1, d), c);
        const auto x3 = linear_sum(scale(x2, d), c);
        const mpq_class want =
            ((vx * exact(d) + vc) * exact(d) + vc) * exact(d) + vc;
        const std::string ex =
            describe(x0) + ", " + describe(c) + ", d = " + show(d);
        exact_ok.add(value(x3) == want, ex);
        form.add(well_formed(x1) && well_formed(x2) && well_formed(x3), ex);
        signs.add(sign(x3) == sgn(want), ex);
    }
    exact_ok.report("chain of scale and linear_sum: exact");
    form.report("chain of scale and linear_sum: well formed at every step");
    signs.report("chain of scale and linear_sum: the right sign");
}

// Four inputs with the components Shewchuk's Linear-Expansion-Sum gives:
// three from the fast_sum story, and one where the other exact sums give
// fewer components. Each was computed by a replay of the algorithm in Python
// floats, and agrees with linear_expansion_sum_zeroelim in Shewchuk's
// predicates.c (both run 2026-10-01).
void test_pinned()
{
    // Shewchuk's own example (section 2.4, in four bits: 11110000 + 1111 +
    // 0.1111 added to itself), in 53 bits. Valid, nonoverlapping, touching
    // throughout. fast_sum's result for it overlaps and even repeats a
    // magnitude, +-2^53, which append() refuses -- so fast_sum is not called
    // here. linear_sum gives two components.
    const double t = std::ldexp(1.0, 53) - 1; // 53 ones
    const Expansion<3> e =
        components<3>({1 - std::ldexp(1.0, -53), t, t * std::ldexp(1.0, 53)});
    check(well_formed(e) && !strongly_nonoverlapping(e) &&
              same_components(linear_sum(e, e),
                              components<2>({-0x1p-52, 0x1p+107})),
          "pinned: Shewchuk's overlap example gives [-2^-52, 2^107]");

    // The counterexamples to Theorem 13: both inputs nonadjacent, and
    // fast_sum gives [1, 6, -2^57] and [-1, 6], where 6 touches the
    // component below it without being a power of two.
    const Expansion<3> a1 = components<3>({1.0, -0x3p+52, -0x1p+57});
    const Expansion<1> b1 = components<1>({0x3p+52 + 6});
    check(nonadjacent(a1) &&
              same_components(linear_sum(a1, b1),
                              components<3>({-1.0, 8.0, -0x1p+57})),
          "pinned: Theorem 13's counterexample gives [-1, 8, -2^57]");
    const Expansion<3> a2 = components<3>({-1.0, -4.0, -0x1p+53});
    const Expansion<1> b2 = components<1>({0x1p+53 + 10});
    check(nonadjacent(a2) &&
              same_components(linear_sum(a2, b2), components<2>({1.0, 4.0})),
          "pinned: the smaller counterexample gives [1, 4]");

    // The value is -3, one double, and sum() and compress(fast_sum()) give
    // just that. Linear-Expansion-Sum does not join components: two ties to
    // even leave a 1 behind, then the carry cancels down to -4. No two
    // components share a magnitude, so the merge order is not in question.
    const Expansion<1> a3 = components<1>({0x3p+52});
    const Expansion<2> b3 = components<2>({-(0x1p+52 + 3), -0x1p+53});
    check(same_components(linear_sum(a3, b3), components<2>({1.0, -4.0})) &&
              sum(a3, b3).size() == 1,
          "pinned: [3 * 2^52] + [-(2^52 + 3), -2^53] gives [1, -4]");
}

// The operators and product() sum with linear_sum: the same components, not
// only the same value. Fails until they are switched over.
void test_switch()
{
    Tally plus, minus, prod;
    std::mt19937_64 rng(96);
    for (int i = 0; i < 20000; ++i)
    {
        mpq_class ignored;
        const Expansion<4> a =
            (i % 2) ? dense<4>(rng, ignored) : build<4>(rng, ignored);
        const Expansion<3> b = dense<3>(rng, ignored);
        const std::string ex = describe(a) + ", " + describe(b);
        plus.add(same_components(a + b, linear_sum(a, b)), ex);
        minus.add(same_components(a - b, linear_sum(a, negate(b))), ex);

        // product() as it should read: the partial products scale(a, b[j]),
        // accumulated with linear_sum, each brought back to the capacity.
        Expansion<24> acc;
        for (std::size_t j = 0; j < b.size(); j++)
            acc = fit<24>(linear_sum(acc, scale(a, b[j])));
        prod.add(same_components(product(a, b), acc), ex);
    }
    plus.report("switch: a + b is linear_sum(a, b)");
    minus.report("switch: a - b is linear_sum(a, negate(b))");
    prod.report("switch: product() accumulates with linear_sum");

    const Expansion<3> a1 = components<3>({1.0, -0x3p+52, -0x1p+57});
    const Expansion<1> b1 = components<1>({0x3p+52 + 6});
    check(same_components(a1 + b1, components<3>({-1.0, 8.0, -0x1p+57})),
          "switch: + on Theorem 13's counterexample gives [-1, 8, -2^57]");
}

} // namespace

int main()
{
    test_grown();
    test_made_other_ways();
    test_dense();
    test_edges();
    test_tree();
    test_chain();
    test_pinned();
    test_switch();
    return report_checks();
}
