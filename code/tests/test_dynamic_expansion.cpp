// DynamicExpansion (kernel/dynamic_expansion.h): expansions sized at runtime,
// and compress, against exact rationals and against the fixed Expansion<N>.
//
// The operations are the fixed ones' algorithms on runtime lengths, so they
// are held to the same contracts -- exact, well formed, and fast_sum's strong
// nonoverlap -- and to the same values as the fixed versions on the same
// inputs. What is new is checked on its own:
//
//   storage   past the 32 inline components, onto the heap and back through
//             copies, which must stay independent of each other
//   compress  exact and nonadjacent, never longer than its input, and with
//             the property plain expansions lack: its largest component alone
//             is within an ulp of the value (Shewchuk's Theorem 23). Not
//             checked, because it is not promised: that the result is the
//             shortest possible, or that compressing again changes nothing --
//             it can. And not visible to these checks: compress's second
//             pass, which a first pass alone already satisfies them without;
//             read it against the paper.
//   operator* compresses, so a formula's intermediate values stay as short as
//             their true values allow
#include <cmath>
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

// The largest component is within 2^-52 of the value, relatively.
bool top_is_close(const D &e, const mpq_class &v)
{
    if (e.size() == 0)
        return v == 0;
    const mpq_class err = abs(exact(e[e.size() - 1]) - v);
    return err <= abs(v) * mpq_class(1, mpz_class(1) << 52);
}

void test_construction_and_storage()
{
    check(D().size() == 0 && sign(D()) == 0, "a default DynamicExpansion is zero");
    check(D(0.0).size() == 0 && D(-3.5).size() == 1 && D(-3.5)[0] == -3.5,
          "from a double: one component, none for zero");
    const TwoTerm t = two_sum(1.0, std::ldexp(1.0, -70));
    check(D(t).size() == 2 && D(t)[0] == t.lo && D(t)[1] == t.hi,
          "from a two-term: the smaller part first");

    std::mt19937_64 rng(80);
    mpq_class v;
    const Expansion<6> fixed = build<6>(rng, v);
    check(value(dyn(fixed)) == v && dyn(fixed).size() == fixed.size(),
          "from an Expansion<N>: the same components");

    // 100 components, one bit each and far apart: well past the inline
    // buffer, so the storage must move to the heap and keep every component.
    D big;
    mpq_class bigv = 0;
    for (int k = 0; k < 100; ++k)
    {
        const double c = std::ldexp(1.0, -500 + 10 * k);
        big.append(c);
        bigv += exact(c);
    }
    check(big.size() == 100 && value(big) == bigv && well_formed(big),
          "100 components: past the inline buffer, every one kept");

    // Copies are independent, inline or on the heap.
    D copy = big;
    copy.append(std::ldexp(1.0, 600));
    check(big.size() == 100 && copy.size() == 101 && value(big) == bigv,
          "appending to a heap-stored copy leaves the original alone");
    D small(2.0);
    D small_copy = small;
    small_copy.append(16.0);
    check(small.size() == 1 && small_copy.size() == 2,
          "appending to an inline copy leaves the original alone");
    D moved = std::move(copy);
    check(moved.size() == 101 && value(moved) == bigv + exact(std::ldexp(1.0, 600)),
          "a moved expansion keeps its components");
}

void test_operations()
{
    Tally exact_ok, form_ok, same_as_fixed;
    std::mt19937_64 rng(81);
    for (int i = 0; i < 20000; ++i)
    {
        mpq_class va, vb;
        const Expansion<6> fa = build<6>(rng, va);
        const Expansion<4> fb = build<4>(rng, vb);
        const D a = dyn(fa), b = dyn(fb);
        const double d = random_double(rng, -100, 100);
        const std::string ex = describe(a) + ", " + describe(b);

        const D g = grow(a, d);
        const D s = fast_sum(a, b);
        const D p = scale(a, d);
        const D q = product(a, b);
        const D n = negate(a);
        exact_ok.add(value(g) == va + exact(d) && value(s) == va + vb &&
                         value(p) == va * exact(d) && value(q) == va * vb &&
                         value(n) == -va,
                     ex);
        form_ok.add(nonadjacent(g) && strongly_nonoverlapping(s) &&
                        nonadjacent(p) && strongly_nonoverlapping(q) &&
                        nonadjacent(n),
                    ex);
        same_as_fixed.add(value(s) == value(fast_sum(fa, fb)) &&
                              value(q) == value(product(fa, fb)) &&
                              sign(q) == sign(product(fa, fb)) &&
                              approximate(s) == approximate(fast_sum(fa, fb)),
                          ex);
    }
    exact_ok.report("dynamic: grow, fast_sum, scale, product, negate are exact");
    form_ok.report("dynamic: each result has its operation's form");
    same_as_fixed.report("dynamic: the same values as the fixed Expansion<N>");

    // Long operands, so that products and sums run past the inline buffer.
    mpq_class va, vb;
    D a = dyn(build<12>(rng, va)), b = dyn(build<12>(rng, vb));
    for (int k = 0; k < 3; ++k)
    {
        a = fast_sum(a, scale(a, 3.0 + k));
        va += va * (3 + k);
    }
    const D q = product(a, b);
    check(value(q) == va * vb && strongly_nonoverlapping(q),
          "a product of long expansions, on the heap, is exact");
}

void test_compress()
{
    Tally exact_ok, form_ok, shorter, close;
    std::mt19937_64 rng(82);
    for (int i = 0; i < 20000; ++i)
    {
        // Long, uncompressed inputs: products of grown expansions.
        mpq_class va, vb;
        const D a = dyn(build<6>(rng, va)), b = dyn(build<5>(rng, vb));
        const D p = product(a, b);
        const mpq_class v = va * vb;
        const D c = compress(p);
        const std::string ex = describe(p);
        exact_ok.add(value(c) == v, ex);
        form_ok.add(nonadjacent(c), ex);
        shorter.add(c.size() <= p.size(), ex);
        close.add(top_is_close(c, v), ex);
    }
    exact_ok.report("compress: the same value");
    form_ok.report("compress: the result is nonadjacent");
    shorter.report("compress: never longer than its input");
    close.report("compress: the largest component is within an ulp of the value");

    check(compress(D()).size() == 0, "compress of zero is zero");
    check(compress(D(2.5)).size() == 1, "compress of one component is itself");
    // Like the example in expansion.h: [1023.5, 1024] is a valid expansion
    // whose largest component is half its value -- but whose value, 2047.5,
    // is a double. Compressed, it is that one double.
    D two;
    two.append(1023.5);
    two.append(1024.0);
    check(compress(two).size() == 1 && compress(two)[0] == 2047.5,
          "compress joins components whose sum is a double");
}

void test_operators()
{
    Tally exact_ok, compressed;
    std::mt19937_64 rng(83);
    for (int i = 0; i < 20000; ++i)
    {
        mpq_class va, vb;
        const D a = dyn(build<4>(rng, va)), b = dyn(build<3>(rng, vb));
        const double d = random_double(rng, -100, 100);
        const std::string ex = describe(a) + ", " + describe(b);
        exact_ok.add(value(a + b) == va + vb && value(a - b) == va - vb &&
                         value(-a) == -va && value(a * b) == va * vb &&
                         value(a * d) == va * exact(d) &&
                         value((a - b) * (a + b) - a * a + b * b) == 0,
                     ex);
        // * compresses: no longer than compress would make the product.
        compressed.add((a * b).size() == compress(product(a, b)).size() &&
                           nonadjacent(a * b),
                       ex);
    }
    exact_ok.report("dynamic operators: +, -, unary -, * are exact");
    compressed.report("dynamic operators: * returns the compressed product");
}

} // namespace

int main()
{
    test_construction_and_storage();
    test_operations();
    test_compress();
    test_operators();
    return report_checks();
}
