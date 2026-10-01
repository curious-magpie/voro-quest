// Expansions (kernel/expansion.h), against exact rationals.
//
// An expansion is a sum of doubles that represents a number exactly. Every
// operation is checked for the two things that make it one:
//
//   exact        its components add up, as rationals, to the true result
//   well formed  its components are nonzero, nonoverlapping, and in order of
//                increasing magnitude -- the invariant every later operation
//                assumes of its inputs, and the reason sign() can look at the
//                largest component alone
//
// and sign() and approximate() are checked against the exact value.
//
// Random inputs are built with grow(), whose own exactness and form are tested
// first, so a fault in it shows up there rather than as a puzzle further down.
#include <cmath>
#include <random>
#include <string>
#include <type_traits>

#include "expansion_support.h"

namespace
{

// --- the operations ---

void test_construction()
{
    const Expansion<3> zero;
    check(zero.size() == 0 && sign(zero) == 0 && approximate(zero) == 0.0,
          "a default expansion is zero, with no components");
    check(Expansion<1>(0.0).size() == 0, "zero as a double has no components");
    check(Expansion<1>(-2.5).size() == 1 && Expansion<1>(-2.5)[0] == -2.5,
          "a nonzero double is one component");

    const TwoTerm t = two_sum(1.0, std::ldexp(1.0, -60));
    const Expansion<2> e(t);
    check(e.size() == 2 && e[0] == t.lo && e[1] == t.hi,
          "a two-term is stored smaller part first");
    check(Expansion<2>(two_sum(1.0, 2.0)).size() == 1,
          "a two-term with no error is one component");

    const Expansion<6> wide(e);
    check(wide.size() == 2 && wide[0] == e[0] && wide[1] == e[1],
          "widening keeps the components");

    // The capacities are part of the type, and follow the operations' bounds.
    static_assert(Expansion<4>::capacity == 4, "capacity is the template size");
    static_assert(std::is_same_v<decltype(grow(Expansion<3>(), 1.0)),
                                 Expansion<4>>,
                  "grow adds one component");
    static_assert(std::is_same_v<decltype(sum(Expansion<3>(), Expansion<5>())),
                                 Expansion<8>>,
                  "sum adds the capacities");
    static_assert(std::is_same_v<decltype(scale(Expansion<3>(), 2.0)),
                                 Expansion<6>>,
                  "scale doubles the capacity");
    static_assert(std::is_same_v<decltype(negate(Expansion<3>())),
                                 Expansion<3>>,
                  "negate keeps the capacity");
}

void test_grow()
{
    Tally exact_ok, form_ok, apart_ok;
    std::mt19937_64 rng(10);
    for (int i = 0; i < 20000; ++i)
    {
        mpq_class want;
        const Expansion<12> e = build<12>(rng, want);
        exact_ok.add(value(e) == want, describe(e));
        form_ok.add(well_formed(e), describe(e));
        apart_ok.add(nonadjacent(e), describe(e));
    }
    exact_ok.report("grow: the components add up to the sum of the inputs");
    form_ok.report("grow: the result is nonzero, nonoverlapping, increasing");
    apart_ok.report("grow: the result is nonadjacent");

    // Exact cancellation, all the way to nothing.
    const Expansion<2> ab = grow(Expansion<1>(1e10), 3e-10);
    const Expansion<4> none = grow(grow(ab, -1e10), -3e-10);
    check(none.size() == 0 && sign(none) == 0,
          "grow cancels exactly to an empty expansion");
}

void test_sum_negate_sign()
{
    Tally exact_sum, form_sum, apart_sum, exact_neg, form_neg, signs, approx;
    std::mt19937_64 rng(11);
    for (int i = 0; i < 20000; ++i)
    {
        mpq_class va, vb;
        const Expansion<7> a = build<7>(rng, va);
        const Expansion<5> b = build<5>(rng, vb);
        const std::string ex = describe(a) + " + " + describe(b);

        const Expansion<12> s = sum(a, b);
        exact_sum.add(value(s) == va + vb, ex);
        form_sum.add(well_formed(s), ex);
        apart_sum.add(nonadjacent(s), ex);

        const Expansion<7> n = negate(a);
        exact_neg.add(value(n) == -va, describe(a));
        form_neg.add(well_formed(n), describe(a));

        // sign and approximate, on the values made above.
        signs.add(sign(s) == sgn(value(s)) && sign(n) == -sgn(va) &&
                      sign(a) == sgn(va),
                  ex);
        const mpq_class vs = va + vb;
        const mpq_class err = abs(mpq_class(approximate(s)) - vs);
        approx.add(err <= abs(vs) * mpq_class(1, mpz_class(1) << 51), ex);
    }
    exact_sum.report("sum: the components add up to a + b");
    form_sum.report("sum: the result is nonzero, nonoverlapping, increasing");
    apart_sum.report("sum: the result is nonadjacent");
    exact_neg.report("negate: the components add up to -a");
    form_neg.report("negate: the result is well formed");
    signs.report("sign: the sign of the exact value");
    approx.report("approximate: within 2^-51 of the exact value, relatively");

    // a + (-a) is exactly zero, with nothing left over.
    mpq_class va;
    const Expansion<7> a = build<7>(rng, va);
    const Expansion<14> z = sum(a, negate(a));
    check(z.size() == 0 && sign(z) == 0, "sum of a and -a is exactly zero");
    check(value(sum(a, Expansion<3>())) == va, "adding zero changes nothing");
}

void test_scale()
{
    Tally exact_ok, form_ok, apart_ok, signs;
    std::mt19937_64 rng(12);
    for (int i = 0; i < 20000; ++i)
    {
        mpq_class va;
        const Expansion<6> a = build<6>(rng, va);
        const double b = random_double(rng, -100, 100);
        const std::string ex = describe(a) + " * " + show(b);

        const Expansion<12> p = scale(a, b);
        exact_ok.add(value(p) == va * exact(b), ex);
        form_ok.add(well_formed(p), ex);
        apart_ok.add(nonadjacent(p), ex);
        signs.add(sign(p) == sgn(va) * (b > 0 ? 1 : -1), ex);
    }
    exact_ok.report("scale: the components add up to a * b");
    form_ok.report("scale: the result is nonzero, nonoverlapping, increasing");
    apart_ok.report("scale: the result is nonadjacent");
    signs.report("scale: the sign is the product of the signs");

    mpq_class va;
    const Expansion<6> a = build<6>(rng, va);
    check(scale(a, 0.0).size() == 0, "scaling by zero gives zero");
    check(value(scale(a, -1.0)) == -va, "scaling by -1 negates");
    check(value(scale(a, 0.5)) == va / 2 && well_formed(scale(a, 0.5)),
          "scaling by a power of two is exact and well formed");
}

} // namespace

int main()
{
    test_construction();
    test_grow();
    test_sum_negate_sign();
    test_scale();
    return report_checks();
}
