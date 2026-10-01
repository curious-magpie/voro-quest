// Expansion arithmetic for formulas (kernel/expansion.h): fit, the product of
// two expansions, and the operators that let one formula be written once and
// evaluated with doubles, expansions or rationals.
//
// Each is checked against exact rationals, for its value and for its form:
// the product goes through fast_sum, so its result is held to fast_sum's
// contract, strongly nonoverlapping.
#include <cmath>
#include <random>
#include <string>
#include <type_traits>

#include "expansion_support.h"

namespace
{

void test_fit()
{
    std::mt19937_64 rng(30);
    mpq_class v;
    const Expansion<6> e = build<6>(rng, v);

    // Narrowing to a capacity that still holds every component.
    const Expansion<6> same = fit<6>(e);
    check(value(same) == v && same.size() == e.size(),
          "fit to the same capacity keeps the value");

    const Expansion<12> wide(e);
    const Expansion<6> back = fit<6>(wide);
    check(value(back) == v, "fit back down from a wider capacity keeps the value");

    static_assert(std::is_same_v<decltype(fit<3>(Expansion<9>())), Expansion<3>>,
                  "fit's result has the capacity asked for");
    check(fit<1>(Expansion<9>()).size() == 0, "fit of zero is zero");
}

void test_product()
{
    Tally exact_ok, form_ok, apart_ok;
    std::mt19937_64 rng(31);
    for (int i = 0; i < 20000; ++i)
    {
        // Within kernel/eft.h's product range: build's inputs are within
        // 2^+-200, and the components of two such expansions multiply to well
        // above 2^-969.
        mpq_class va, vb;
        const Expansion<4> a = build<4>(rng, va);
        const Expansion<3> b = build<3>(rng, vb);
        const Expansion<24> p = product(a, b);
        const std::string ex = describe(a) + " * " + describe(b);
        exact_ok.add(value(p) == va * vb, ex);
        form_ok.add(well_formed(p), ex);
        apart_ok.add(strongly_nonoverlapping(p), ex);
    }
    exact_ok.report("product: exactly a * b");
    form_ok.report("product: well formed");
    apart_ok.report("product: strongly nonoverlapping");

    static_assert(std::is_same_v<decltype(product(Expansion<4>(), Expansion<3>())),
                                 Expansion<24>>,
                  "product has capacity 2 * M * N");

    mpq_class va;
    const Expansion<4> a = build<4>(rng, va);
    check(product(a, Expansion<2>()).size() == 0, "a product with zero is zero");
    check(value(product(a, Expansion<1>(-1.0))) == -va, "a product with -1 negates");
    check(value(product(Expansion<1>(3.0), Expansion<1>(0.1))) ==
              exact(3.0) * exact(0.1),
          "the product of two doubles is exact");
}

void test_operators()
{
    Tally exact_ok;
    std::mt19937_64 rng(32);
    for (int i = 0; i < 20000; ++i)
    {
        mpq_class va, vb;
        const Expansion<3> a = build<3>(rng, va);
        const Expansion<2> b = build<2>(rng, vb);
        const double d = random_double(rng, -100, 100);
        const std::string ex = describe(a) + ", " + describe(b);

        exact_ok.add(value(a + b) == va + vb, ex + " +");
        exact_ok.add(value(a - b) == va - vb, ex + " -");
        exact_ok.add(value(-a) == -va, ex + " unary -");
        exact_ok.add(value(a * b) == va * vb, ex + " *");
        exact_ok.add(value(a * d) == va * exact(d), ex + " * double");
        // A small formula, written the way a predicate is.
        exact_ok.add(value((a - b) * (a + b) - a * a + b * b) == 0, ex + " formula");
    }
    exact_ok.report("operators: +, -, unary -, * are exact");

    static_assert(std::is_same_v<decltype(Expansion<3>() + Expansion<2>()),
                                 Expansion<5>>,
                  "+ is fast_sum");
    static_assert(std::is_same_v<decltype(Expansion<3>() - Expansion<2>()),
                                 Expansion<5>>,
                  "- is fast_sum with a negation");
    static_assert(std::is_same_v<decltype(Expansion<3>() * Expansion<2>()),
                                 Expansion<12>>,
                  "* of two expansions is product");
    static_assert(std::is_same_v<decltype(Expansion<3>() * 2.0), Expansion<6>>,
                  "* by a double is scale");
}

} // namespace

int main()
{
    test_fit();
    test_product();
    test_operators();
    return report_checks();
}
