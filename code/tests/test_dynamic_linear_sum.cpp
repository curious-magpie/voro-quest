// linear_sum for DynamicExpansion (kernel/dynamic_expansion.h): the same
// Linear-Expansion-Sum as expansion.h's (test_linear_sum.cpp), on runtime
// lengths. Write that one first: this suite checks against it.
//
// The contract is the same, Theorem 24 -- nonoverlapping inputs, an exact
// nonoverlapping result -- and so are the checks, plus what is new here:
//
//   fixed      the same components as the Expansion<N> version, input for
//              input: one algorithm, two storages
//   long       inputs and results past the 32 inline components
//   switch     a + b, a - b and product() are linear_sum, and * compresses
//              only its result: the compress of the left factor, there to
//              make fast_sum's inputs nonadjacent, has no job left, so
//              a * b is compress(product(a, b)) and a * d is scale(a, d).
//              These checks fail until dynamic_expansion.h is switched over
#include <cmath>
#include <initializer_list>
#include <random>
#include <string>

#include "expansion_support.h"
#include "kernel/dynamic_expansion.h"

namespace
{

using D = DynamicExpansion;

template <std::size_t N> D dyn(const Expansion<N> &e)
{
    return D(e);
}

D components(std::initializer_list<double> cs)
{
    D e;
    for (double c : cs)
        e.append(c);
    return e;
}

struct Checks
{
    Tally exact, form, signs;

    void run(const D &a, const D &b)
    {
        const D s = linear_sum(a, b);
        const std::string ex = describe(a) + " + " + describe(b);
        const mpq_class want = value(a) + value(b);
        exact.add(value(s) == want, ex);
        form.add(well_formed(s), ex);
        signs.add(sign(s) == sgn(want), ex);
    }

    void report(const std::string &what) const
    {
        exact.report(
            ("dynamic linear_sum, " + what + ": exactly a + b").c_str());
        form.report(("dynamic linear_sum, " + what + ": well formed").c_str());
        signs.report(
            ("dynamic linear_sum, " + what + ": the right sign").c_str());
    }
};

void test_same_as_fixed()
{
    Checks c;
    Tally same;
    std::mt19937_64 rng(100);
    mpq_class ignored;
    for (int i = 0; i < 20000; ++i)
    {
        const Expansion<6> fa =
            (i % 2) ? dense<6>(rng, ignored) : build<6>(rng, ignored);
        const Expansion<4> fb = dense<4>(rng, ignored);
        const Expansion<1> f1 = build<1>(rng, ignored);
        const D a = dyn(fa), b = dyn(fb), one = dyn(f1);
        c.run(a, b);
        c.run(b, a);
        c.run(a, a);
        c.run(a, negate(a));
        c.run(a, one);
        c.run(one, one);
        const std::string ex = describe(a) + " + " + describe(b);
        same.add(same_components(linear_sum(a, b), linear_sum(fa, fb)) &&
                     same_components(linear_sum(b, a), linear_sum(fb, fa)) &&
                     same_components(linear_sum(a, one), linear_sum(fa, f1)) &&
                     same_components(linear_sum(one, one), linear_sum(f1, f1)),
                 ex);
    }
    c.report("inputs of every shape");
    same.report("dynamic linear_sum: the same components as the fixed one");
}

void test_edges()
{
    Checks c;
    std::mt19937_64 rng(101);
    mpq_class ignored;
    const D a = dyn(dense<6>(rng, ignored));
    const D zero;
    c.run(zero, zero);
    c.run(zero, a);
    c.run(a, zero);
    c.run(zero, D(3.0));
    c.run(D(3.0), zero);
    c.run(D(1.0), D(-1.0));
    c.run(D(std::ldexp(1.0, 100)), D(1.0));
    c.report("edge cases");

    check(linear_sum(zero, zero).size() == 0,
          "dynamic linear_sum of zeros has no components");
    check(linear_sum(a, negate(a)).size() == 0,
          "dynamic linear_sum of a and -a has no components");
    check(same_components(linear_sum(zero, D(3.0)), D(3.0)) &&
              same_components(linear_sum(D(3.0), zero), D(3.0)),
          "dynamic linear_sum of one component and zero is that component");
}

// Past the inline buffer: 40 and 30 components in, up to 70 out.
void test_long()
{
    Checks c;
    std::mt19937_64 rng(102);
    for (int i = 0; i < 2000; ++i)
    {
        mpq_class va, vb;
        const D a = dyn(dense<40>(rng, va, 10));
        const D b = dyn(dense<30>(rng, vb, 10));
        c.run(a, b);
        c.run(a, a);
        c.run(b, a);
    }
    c.report("long inputs");
}

// x -> x * d + c, six times, at runtime lengths: each step's input is the
// last step's output, and each step is a theorem.
void test_chain()
{
    Tally exact_ok, form, signs;
    std::mt19937_64 rng(103);
    for (int i = 0; i < 10000; ++i)
    {
        mpq_class vx, vc;
        D x = dyn(dense<4>(rng, vx, 20));
        const D c = dyn(dense<3>(rng, vc, 20));
        const double d = (i % 2) ? random_double(rng, -20, 20)
                                 : double(int(rng() % 31) - 15);
        mpq_class want = vx;
        bool formed = true;
        for (int k = 0; k < 6; ++k)
        {
            x = linear_sum(scale(x, d), c);
            want = want * exact(d) + vc;
            formed = formed && well_formed(x);
        }
        const std::string ex = describe(c) + ", d = " + show(d);
        exact_ok.add(value(x) == want, ex);
        form.add(formed, ex);
        signs.add(sign(x) == sgn(want), ex);
    }
    exact_ok.report("dynamic chain of scale and linear_sum: exact");
    form.report("dynamic chain of scale and linear_sum: well formed");
    signs.report("dynamic chain of scale and linear_sum: the right sign");
}

void test_pinned()
{
    const double t = std::ldexp(1.0, 53) - 1;
    const D e =
        components({1 - std::ldexp(1.0, -53), t, t * std::ldexp(1.0, 53)});
    check(same_components(linear_sum(e, e), components({-0x1p-52, 0x1p+107})),
          "dynamic pinned: Shewchuk's overlap example gives [-2^-52, 2^107]");
    const D a1 = components({1.0, -0x3p+52, -0x1p+57});
    const D b1 = components({0x3p+52 + 6});
    check(
        same_components(linear_sum(a1, b1), components({-1.0, 8.0, -0x1p+57})),
        "dynamic pinned: Theorem 13's counterexample gives [-1, 8, -2^57]");
    const D a2 = components({-1.0, -4.0, -0x1p+53});
    const D b2 = components({0x1p+53 + 10});
    check(same_components(linear_sum(a2, b2), components({1.0, 4.0})),
          "dynamic pinned: the smaller counterexample gives [1, 4]");
    // One double's worth (-3) in two components; compress joins them.
    const D a3 = components({0x3p+52});
    const D b3 = components({-(0x1p+52 + 3), -0x1p+53});
    check(same_components(linear_sum(a3, b3), components({1.0, -4.0})) &&
              same_components(compress(linear_sum(a3, b3)), D(-3.0)),
          "dynamic pinned: [3 * 2^52] + [-(2^52 + 3), -2^53] gives [1, -4]");
}

void test_switch()
{
    Tally plus, minus, prod, times, scaled;
    std::mt19937_64 rng(104);
    for (int i = 0; i < 20000; ++i)
    {
        mpq_class va, vb;
        const D a = (i % 2) ? dyn(dense<5>(rng, va)) : dyn(build<5>(rng, va));
        const D b = dyn(dense<4>(rng, vb));
        const std::string ex = describe(a) + ", " + describe(b);
        plus.add(same_components(a + b, linear_sum(a, b)), ex);
        minus.add(same_components(a - b, linear_sum(a, negate(b))), ex);
        // product() as it should read: no fit, no compress.
        D acc;
        for (std::size_t j = 0; j < b.size(); j++)
            acc = linear_sum(acc, scale(a, b[j]));
        prod.add(same_components(product(a, b), acc), ex);
        // * compresses the product, and nothing else.
        times.add(same_components(a * b, compress(product(a, b))) &&
                      value(a * b) == va * vb,
                  ex);
        const double d = random_double(rng, -100, 100);
        scaled.add(same_components(a * d, scale(a, d)), ex);
    }
    plus.report("dynamic switch: a + b is linear_sum(a, b)");
    minus.report("dynamic switch: a - b is linear_sum(a, negate(b))");
    prod.report("dynamic switch: product() accumulates with linear_sum");
    times.report("dynamic switch: a * b is compress(product(a, b))");
    scaled.report("dynamic switch: a * d is scale(a, d)");
}

} // namespace

int main()
{
    test_same_as_fixed();
    test_edges();
    test_long();
    test_chain();
    test_pinned();
    test_switch();
    return report_checks();
}
