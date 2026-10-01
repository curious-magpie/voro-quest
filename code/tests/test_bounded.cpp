// Bounded (kernel/bounded.h): a double with a certified bound on its error, the
// number type of the filters. Checked against exact rationals.
//
// The one promise that matters is soundness: the true value of every formula
// lies within error of value, so a sign the filter calls certain is the true
// sign. A filter that is not sound gives wrong answers silently, which is the
// worst thing a predicate can do, so it is checked on formulas of every shape
// the predicates use, with inputs chosen to cancel.
//
// The second promise is usefulness: the bound is not much larger than it has
// to be, or the filter would hand everything to the exact path. That is
// checked on single operations, where the rounding is known.
#include <cmath>
#include <optional>
#include <random>
#include <string>

#include "test_support.h"

#include "exact.h"
#include "kernel/bounded.h"

namespace
{

using B = Bounded;

// |value - exact| <= error, as rationals, and the error a finite non-negative
// double.
bool sound(const B &b, const mpq_class &want)
{
    if (!(b.error >= 0.0) || !std::isfinite(b.error) || !std::isfinite(b.value))
        return false;
    return abs(exact(b.value) - want) <= exact(b.error);
}

// A certain sign must be the true one; no sign is always allowed.
bool sign_honest(const B &b, const mpq_class &want)
{
    const std::optional<int> s = certain_sign(b);
    return !s || *s == sgn(want);
}

struct Tally
{
    long failures = 0;
    std::string first;

    void add(bool ok, const std::string &example)
    {
        if (!ok && failures++ == 0)
            first = example;
    }
    void report(const char *what) const
    {
        check(failures == 0, what, failures ? first.c_str() : nullptr);
    }
};

std::string describe(const B &b, const mpq_class &want)
{
    return "value " + show(b.value) + ", error " + show(b.error) + ", exact " +
           show(want.get_d());
}

// Eight doubles near one another, so that differences cancel to a few bits,
// or spread out, so that they do not.
void draw(std::mt19937_64 &rng, double (&x)[8], bool close)
{
    const double base = random_double(rng, -30, 30);
    for (double &v : x)
        v = close ? base * (1.0 + random_double(rng, -50, -20))
                  : random_double(rng, -30, 30);
}

void test_soundness()
{
    Tally ok, honest;
    std::mt19937_64 rng(50);
    for (int i = 0; i < 100000; ++i)
    {
        double x[8];
        draw(rng, x, i % 2 == 0);
        B b[8];
        mpq_class q[8];
        for (int k = 0; k < 8; ++k)
        {
            b[k] = B(x[k]);
            q[k] = exact(x[k]);
        }

        // The shape of orient2d: two products of differences, subtracted.
        const B d = (b[0] - b[1]) * (b[2] - b[3]) - (b[4] - b[5]) * (b[6] - b[7]);
        const mpq_class dq =
            (q[0] - q[1]) * (q[2] - q[3]) - (q[4] - q[5]) * (q[6] - q[7]);
        ok.add(sound(d, dq), describe(d, dq));
        honest.add(sign_honest(d, dq), describe(d, dq));

        // A chain of products, where relative errors pile up.
        B p = b[0];
        mpq_class pq = q[0];
        for (int k = 1; k < 8; ++k)
        {
            p = p * b[k];
            pq *= q[k];
        }
        ok.add(sound(p, pq), describe(p, pq));
        honest.add(sign_honest(p, pq), describe(p, pq));

        // A chain of sums that cancels, then is multiplied: errors from the
        // cancelled part are all that is left.
        B s = b[0] + b[1] - b[0] - b[1] + b[2];
        const mpq_class sq = q[2];
        s = -(s * (b[3] - b[4]));
        const mpq_class tq = -(sq * (q[3] - q[4]));
        ok.add(sound(s, tq), describe(s, tq));
        honest.add(sign_honest(s, tq), describe(s, tq));

        // A deeper mix: the shape of a 3 x 3 determinant's cofactor terms.
        const B m = ((b[0] + b[1]) * (b[2] - b[3]) + b[4]) * (b[5] - b[6]) -
                    b[7] * b[7];
        const mpq_class mq =
            ((q[0] + q[1]) * (q[2] - q[3]) + q[4]) * (q[5] - q[6]) - q[7] * q[7];
        ok.add(sound(m, mq), describe(m, mq));
        honest.add(sign_honest(m, mq), describe(m, mq));
    }
    ok.report("Bounded: the exact value is always within error of value");
    honest.report("Bounded: a certain sign is always the true sign");
}

void test_exact_cases()
{
    // Inputs are exact, so their error is zero, and so is that of anything
    // built only from exact zeros.
    const B a(2.5);
    check(a.value == 2.5 && a.error == 0.0, "an input is exact: error 0");
    check(certain_sign(a) == std::optional<int>(1) &&
              certain_sign(B(-1e-300)) == std::optional<int>(-1),
          "an input's sign is certain");

    const B z = B(0.1) - B(0.1);
    check(z.value == 0.0 && z.error == 0.0 && certain_sign(z) == std::optional<int>(0),
          "x - x is exactly zero, and certainly so");
    check(certain_sign(B(0.0) * B(7.0)) == std::optional<int>(0),
          "a product with an exact zero is certainly zero");
    check((-a).value == -2.5 && (-a).error == 0.0, "negation is exact");
}

void test_tightness()
{
    // One operation on exact inputs rounds once: its error bound should be
    // about half an ulp, and certainly within one, so that the filter is not
    // giving up decisions it could make. 2^-52 relative is one ulp at most.
    Tally sum_tight, product_tight;
    std::mt19937_64 rng(51);
    const double one_ulp = std::ldexp(1.0, -52);
    for (int i = 0; i < 100000; ++i)
    {
        const double x = random_double(rng, -100, 100);
        const double y = random_double(rng, -100, 100);
        const B s = B(x) + B(y), p = B(x) * B(y);
        sum_tight.add(s.error <= one_ulp * std::fabs(s.value), describe(s, 0));
        product_tight.add(p.error <= one_ulp * std::fabs(p.value),
                          describe(p, 0));
    }
    sum_tight.report("Bounded: a sum of exact inputs has an error within an ulp");
    product_tight.report(
        "Bounded: a product of exact inputs has an error within an ulp");
}

} // namespace

int main()
{
    test_soundness();
    test_exact_cases();
    test_tightness();
    return report_checks();
}
