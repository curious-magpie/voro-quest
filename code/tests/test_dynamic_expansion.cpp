// DynamicExpansion (kernel/dynamic_expansion.h): expansions sized at runtime,
// and compress, against exact rationals and against the fixed Expansion<N>.
//
// The operations are the fixed ones' algorithms on runtime lengths, so they
// are held to the same contracts -- exact, well formed, nonadjacent where
// Shewchuk proves it -- and to the same values as the fixed versions on the
// same inputs. What is new is checked on its own:
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
    check(D().size() == 0 && sign(D()) == 0,
          "a default DynamicExpansion is zero");
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
    check(moved.size() == 101 &&
              value(moved) == bigv + exact(std::ldexp(1.0, 600)),
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
        form_ok.add(nonadjacent(g) && well_formed(s) && nonadjacent(p) &&
                        well_formed(q) && nonadjacent(n),
                    ex);
        same_as_fixed.add(value(s) == value(fast_sum(fa, fb)) &&
                              value(q) == value(product(fa, fb)) &&
                              sign(q) == sign(product(fa, fb)) &&
                              approximate(s) == approximate(fast_sum(fa, fb)),
                          ex);
    }
    exact_ok.report(
        "dynamic: grow, fast_sum, scale, product, negate are exact");
    form_ok.report("dynamic: each result has its operation's form");
    same_as_fixed.report("dynamic: the same values as the fixed Expansion<N>");

    // Long operands, so that products and sums run past the inline buffer.
    mpq_class va, vb;
    D a = dyn(build<12>(rng, va)), b = dyn(build<12>(rng, vb));
    for (int k = 0; k < 3; ++k)
    {
        a = a + scale(a, 3.0 + k);
        va += va * (3 + k);
    }
    const D q = product(a, b);
    check(value(q) == va * vb && well_formed(q),
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
    close.report(
        "compress: the largest component is within an ulp of the value");

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

// +-2^e, the sign at random.
double signed_power(std::mt19937_64 &rng, int e)
{
    return (rng() & 1) ? -std::ldexp(1.0, e) : std::ldexp(1.0, e);
}

// A strongly nonoverlapping expansion that is *not* nonadjacent: at least one
// pair of touching powers of two, +-2^k and +-2^(k+1), the one kind of
// adjacency strong nonoverlap allows. The other components have up to 8
// bits, and every gap but a pair's is at least one bit, so no component
// touches two neighbours. This is what a fast_sum result may look like, and
// what random inputs almost never produce: a grown expansion is nonadjacent.
D touching_pairs(std::mt19937_64 &rng, mpq_class &sum)
{
    D e;
    int low = -100 + int(rng() % 20); // lowest bit of the next component
    const int n = 2 + int(rng() % 5);
    bool paired = false;
    for (int i = 0; i < n; ++i)
    {
        if ((rng() & 1) || (i == n - 1 && !paired))
        {
            e.append(signed_power(rng, low));
            e.append(signed_power(rng, low + 1));
            low += 2;
            paired = true;
        }
        else
        {
            // An odd k-bit integer, so its bits are exactly low .. low+k-1.
            const int k = 1 + int(rng() % 8);
            const double m =
                double((1u << (k - 1)) | (rng() & ((1u << k) - 1)) | 1u);
            e.append((rng() & 1) ? -std::ldexp(m, low) : std::ldexp(m, low));
            low += k;
        }
        low += 1 + int(rng() % 4);
    }
    sum = value(e);
    return e;
}

// An expansion from its components, smallest first.
D components(std::initializer_list<double> cs)
{
    D e;
    for (double c : cs)
        e.append(c);
    return e;
}

// The inputs that Shewchuk's theorems leave open: strongly nonoverlapping,
// with touching components, as fast_sum results may be. Theorem 19 promises
// scale() only nonoverlapping output for those, and these inputs show that
// strong nonoverlap can indeed be lost -- by scale(), and by fast_sum itself
// (the two pinned examples at the end). So the checks here are of three
// kinds, labelled as such:
//
//   evidence  what a predicate needs, on chains the theorems do not cover:
//             exact, nonoverlapping in increasing magnitude (well_formed),
//             and so the right sign. Not a theorem; a failure here is a
//             counterexample that matters, and is worth keeping. The chains
//             call fast_sum by name, so they stay a record of fast_sum now
//             that the operators sum with linear_sum.
//   contract  a * b ends in compress, so it is nonadjacent (Theorem 23);
//             a * d is a scale, so it is well formed (Theorem 19) but, for
//             these inputs, not necessarily nonadjacent.
//   pinned    the two counterexamples to strong nonoverlap found by this
//             test (2026-10-01), kept as fixed inputs so they stay on record.
//             If one of these starts failing, the arithmetic changed: look.
void test_touching_inputs()
{
    Tally generator, exact_ok, evidence, operators;
    std::mt19937_64 rng(84);
    for (int i = 0; i < 100000; ++i)
    {
        mpq_class va, vc;
        const D a = touching_pairs(rng, va), c = touching_pairs(rng, vc);
        // Few-bit factors as well as full ones: they make carries line up.
        const double d = (i % 2) ? random_double(rng, -100, 100)
                                 : double(int(rng() % 31) - 15);
        const std::string ex =
            describe(a) + ", " + describe(c) + ", d = " + show(d);
        generator.add(strongly_nonoverlapping(a) && !nonadjacent(a), ex);

        const D s = scale(a, d);
        const D f = fast_sum(s, c);
        const D p = product(a, c);
        // A chain of raw sums of raw scales, never compressed.
        D chain = f;
        mpq_class vchain = va * exact(d) + vc;
        for (int k = 0; k < 4; ++k)
        {
            chain = fast_sum(scale(chain, d), c);
            vchain = vchain * exact(d) + vc;
        }
        exact_ok.add(value(s) == va * exact(d) && value(f) == value(s) + vc &&
                         value(p) == va * vc && value(chain) == vchain &&
                         value(a * c) == va * vc &&
                         value(a * d) == va * exact(d),
                     ex);
        evidence.add(well_formed(s) && well_formed(f) && well_formed(p) &&
                         well_formed(chain) && sign(chain) == sgn(vchain),
                     ex);
        operators.add(nonadjacent(a * c) && well_formed(a * d), ex);
    }
    generator.report("touching inputs: strongly nonoverlapping, not "
                     "nonadjacent");
    exact_ok.report("touching inputs: scale, fast_sum, product, * are exact");
    evidence.report("evidence: raw chains stay nonoverlapping, right sign");
    operators.report("contract: a * b nonadjacent, a * d well formed");

    // scale() of a strongly nonoverlapping expansion (three touching pairs of
    // powers of two) by a 53-bit double: the result has a power of two,
    // 2^-64, touching a 52-bit component. Nonoverlapping, as Theorem 19
    // promises, but not strongly.
    const D a1 =
        components({0x1p-82, -0x1p-81, -0x1p-79, -0x1p-78, -0x1p-76, 0x1p-75});
    const double d1 = 36256385392834236416.0;
    const D s1 = scale(a1, d1);
    check(strongly_nonoverlapping(a1) && well_formed(s1) &&
              value(s1) == value(a1) * exact(d1) &&
              !strongly_nonoverlapping(s1),
          "pinned: scale can lose strong nonoverlap");

    // fast_sum of a nonadjacent and a strongly nonoverlapping expansion --
    // both valid inputs for Theorem 13 -- whose result has -2^-141 adjacent
    // to a 5-bit component (bits -140 .. -136). Traced by hand: the carry
    // drops -2^-141 as an error, then cancels against the 53-bit component
    // down to a lowest bit of -140, and its next error lands just above.
    // Theorem 13 promises a strongly nonoverlapping result. Checked against
    // the paper's statement and section 2.3's definitions (read 2026-10-01),
    // and replayed independently with Python floats: every hypothesis holds
    // and the conclusion does not, so this is a counterexample to Theorem 13
    // as stated. The merge has no ties, so the output is the algorithm's.
    const D s2 = components({8.5528470722950261e-48, -3.941151930913548e-46,
                             -1.4012984643248171e-44, -3.1389085600875902e-43,
                             4.4881243107836986e-27});
    const D c2 = components({-0x1p-89, -0x1p-88, -3.1019272970738538e-25});
    const D f2 = fast_sum(s2, c2);
    check(nonadjacent(s2) && strongly_nonoverlapping(c2) && well_formed(f2) &&
              value(f2) == value(s2) + value(c2) &&
              !strongly_nonoverlapping(f2),
          "pinned: fast_sum can lose strong nonoverlap");
}

} // namespace

int main()
{
    test_construction_and_storage();
    test_operations();
    test_compress();
    test_operators();
    test_touching_inputs();
    return report_checks();
}
