// The error-free transformations of kernel/eft.h, against exact rationals.
//
// Each transformation turns one rounded operation into two doubles whose sum is
// the exact result. Three things are checked for every input:
//
//   exact       hi + lo equals the true a + b (or a * b), as rationals
//   rounded     hi is what the plain double operation returns, bit for bit
//   separated   lo is at most half an ulp of hi, so hi + lo rounds back to hi:
//               the two parts do not overlap, which every expansion relies on
//
// on hand-picked cases where rounding is at its most delicate, and on hundreds
// of thousands of random doubles spread over wide exponent ranges.
#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "test_support.h"

#include "exact.h"
#include "kernel/eft.h"

namespace
{

// Counts the failures of one property over many inputs, and keeps the first,
// so that one check() reports the property with an example to reproduce.
struct Tally
{
    int failures = 0;
    std::string first;

    void add(bool ok, double a, double b)
    {
        if (ok)
            return;
        if (failures++ == 0)
            first = "a = " + show(a) + ", b = " + show(b);
    }
    void report(const char *what) const
    {
        check(failures == 0, what, failures ? first.c_str() : nullptr);
    }
};

// hi + lo rounds back to hi: |lo| is within half an ulp of hi.
bool separated(const TwoTerm &t)
{
    return t.hi + t.lo == t.hi;
}

// Where addition is hard: exact cancellation, ties that round to even, a term
// too small to change the sum, and exponents far apart.
std::vector<std::pair<double, double>> hard_sums()
{
    const double u = std::ldexp(1.0, -53); // half an ulp of 1
    return {
        {1.0, u},                        // a tie: rounds to even, down to 1
        {1.0 + 2 * u, u},                // a tie: rounds to even, up
        {1.0, std::ldexp(1.0, -54)},     // below the tie: lost entirely
        {1.0, -u},                       // just below 1
        {std::ldexp(1.0, 100), 1.0},     // absorbed
        {std::ldexp(1.0, 100), -1.0},    // absorbed, from below
        {3.0, -3.0},                     // exact cancellation to zero
        {0.1, -0.1},                     // the same, off a power of two
        {1.0 + 2 * u, -1.0},             // cancellation that is exact
        {0.0, 1.5},                      // zero on either side
        {-2.5, 0.0},                     //
        {0.0, -0.0},                     // signed zeros
        {1e300, 1e300},                  // large, but not overflowing
        {std::ldexp(1.0, -1000), 1.0},   // tiny against 1
    };
}

std::vector<std::pair<double, double>> hard_products()
{
    const double e = std::ldexp(1.0, -52); // an ulp of 1
    return {
        {1.0 + e, 1.0 + e},              // error term 2^-104
        {1.0 + e, 1.0 - e},              // 1 - 2^-104: the error is all of it
        {3.0, 1.0 / 3.0},                // rounds to 1, error nonzero
        {0.1, 10.0},                     //
        {std::ldexp(1.0, 60), 7.0},      // exact: the error term is zero
        {-(1.0 + e), 1.0 + 2 * e},       // a negative product
        {0.0, 123.0},                    // zero
        {1e150, 1e150},                  // large, but not overflowing
        {1e-140, 1e-140},                // small, but the error term (about
                                         // 1e-296) stays above the subnormals
    };
}

void test_two_sum()
{
    Tally exact_sum, rounded, apart;
    const auto one = [&](double a, double b)
    {
        const TwoTerm t = two_sum(a, b);
        exact_sum.add(exact(t.hi) + exact(t.lo) == exact(a) + exact(b), a, b);
        rounded.add(t.hi == a + b, a, b);
        apart.add(separated(t), a, b);
    };

    for (const auto &[a, b] : hard_sums())
    {
        one(a, b);
        one(b, a); // two_sum has no order precondition, so both must work
    }

    std::mt19937_64 rng(1);
    for (int i = 0; i < 200000; ++i)
    {
        // Far-apart exponents, then close ones, where cancellation happens.
        one(random_double(rng, -500, 500), random_double(rng, -500, 500));
        const double a = random_double(rng, -20, 20);
        one(a, -a * (1.0 + random_double(rng, -60, -20)));
    }

    exact_sum.report("two_sum: hi + lo is exactly a + b");
    rounded.report("two_sum: hi is the rounded sum");
    apart.report("two_sum: lo is within half an ulp of hi");
}

void test_fast_two_sum()
{
    // Its precondition is |a| >= |b|, so every pair is put in that order.
    Tally exact_sum, rounded, apart;
    const auto one = [&](double a, double b)
    {
        if (std::fabs(a) < std::fabs(b))
            std::swap(a, b);
        const TwoTerm t = fast_two_sum(a, b);
        exact_sum.add(exact(t.hi) + exact(t.lo) == exact(a) + exact(b), a, b);
        rounded.add(t.hi == a + b, a, b);
        apart.add(separated(t), a, b);
    };

    for (const auto &[a, b] : hard_sums())
        one(a, b);

    std::mt19937_64 rng(2);
    for (int i = 0; i < 200000; ++i)
    {
        one(random_double(rng, -500, 500), random_double(rng, -500, 500));
        const double a = random_double(rng, -20, 20);
        one(a, -a * (1.0 + random_double(rng, -60, -20)));
    }

    exact_sum.report("fast_two_sum: hi + lo is exactly a + b");
    rounded.report("fast_two_sum: hi is the rounded sum");
    apart.report("fast_two_sum: lo is within half an ulp of hi");
}

void test_two_product()
{
    // Exponents kept within +-250, so that no product overflows and no error
    // term underflows -- the range eft.h promises exactness on.
    Tally exact_fma, rounded_fma, apart_fma, exact_dekker, same;
    const auto one = [&](double a, double b)
    {
        const TwoTerm t = two_product(a, b);
        exact_fma.add(exact(t.hi) + exact(t.lo) == exact(a) * exact(b), a, b);
        rounded_fma.add(t.hi == a * b, a, b);
        apart_fma.add(separated(t), a, b);

        // The rounded product and its error are unique, so the two methods
        // must agree exactly; only a zero error may differ in its sign.
        const TwoTerm d = two_product_dekker(a, b);
        exact_dekker.add(exact(d.hi) + exact(d.lo) == exact(a) * exact(b), a, b);
        same.add(d.hi == t.hi && d.lo == t.lo, a, b);
    };

    for (const auto &[a, b] : hard_products())
    {
        one(a, b);
        one(b, a);
    }

    std::mt19937_64 rng(3);
    for (int i = 0; i < 200000; ++i)
        one(random_double(rng, -250, 250), random_double(rng, -250, 250));

    exact_fma.report("two_product: hi + lo is exactly a * b");
    rounded_fma.report("two_product: hi is the rounded product");
    apart_fma.report("two_product: lo is within half an ulp of hi");
    exact_dekker.report("two_product_dekker: hi + lo is exactly a * b");
    same.report("two_product_dekker agrees with two_product");
}

void test_split()
{
    // Each half must fit in 26 bits, so that a product of two halves is exact
    // in a double -- the whole point of splitting.
    Tally exact_split, narrow;
    const auto one = [&](double a)
    {
        const TwoTerm s = split(a);
        exact_split.add(exact(s.hi) + exact(s.lo) == exact(a), a, 0.0);
        narrow.add(significant_bits(s.hi) <= 26 && significant_bits(s.lo) <= 26,
                   a, 0.0);
    };

    for (double a : {1.0, -1.0, 0.0, 0.1, 1.0 / 3.0, 1.0 + std::ldexp(1.0, -52),
                     std::ldexp(1.0, 900), -std::ldexp(3.0, -900)})
        one(a);

    std::mt19937_64 rng(4);
    for (int i = 0; i < 200000; ++i)
        one(random_double(rng, -900, 900));

    exact_split.report("split: hi + lo is exactly a");
    narrow.report("split: both halves fit in 26 bits");
}

} // namespace

int main()
{
    test_two_sum();
    test_fast_two_sum();
    test_two_product();
    test_split();
    return report_checks();
}
