// Expansions: exact real numbers, as short sums of doubles.
//
// A double holds 53 significant bits. The value of a predicate -- a
// determinant of input coordinates -- can need hundreds, and its sign, which
// is all a predicate wants, can depend on the very last of them. An expansion
// holds as many bits as it needs, as a list of doubles whose sum is the exact
// value:
//
//   x = c[0] + c[1] + ... + c[n-1]        exactly, in the real numbers
//
// built only from the error-free transformations of eft.h, so that nothing is
// ever rounded away: every operation here returns the exact sum, product or
// negation of its operands, as another expansion.
//
// The idea and the algorithms are Shewchuk's (Adaptive precision floating-
// point arithmetic and fast robust geometric predicates, Discrete &
// Computational Geometry 18, 1997, sections 2.3 to 2.6). The theorem numbers
// below are his; check them against the paper when citing.
//
// --- the invariant ---
//
// Every Expansion, whatever made it, has components that are
//
//   nonzero                zero is the empty expansion, never a component
//   nonoverlapping         in the sense of eft.h: each component's lowest set
//                          bit is above the next smaller one's highest
//   increasing magnitude   smallest first: |c[0]| < |c[1]| < ... < |c[n-1]|
//
// Every operation may assume this of its inputs and must keep it in its
// result. Two things follow, and the rest of the kernel relies on both:
//
//   sign is cheap. All the components below the largest add up to less than
//   the largest's lowest set bit, so they cannot change its sign: the sign of
//   the whole is the sign of c[n-1]. A predicate's answer is one comparison.
//
//   But c[n-1] alone is *not* a good approximation of the value. The bound
//   above says the rest is smaller than c[n-1], not that it is small: [1023.9,
//   1024] is nonoverlapping (1023.9 lies wholly below 1024's only set bit), and
//   its value is twice its largest component. To read a value, add the
//   components up: approximate() does, and is accurate.
//
// The invariant is kept in one place, append(), the only way a component gets
// into an Expansion: it drops zeros, and asserts that each component is
// larger than the last. tests/test_expansion.cpp checks the full
// nonoverlapping property, against exact rationals, after every operation.
//
// --- the capacity ---
//
// An Expansion<N> holds at most N components, in a std::array on the stack:
// no allocation, ever, which matters because predicates run millions of
// times. How many components a result can need is known from its operands --
// a + b needs at most as many as a and b together, a * b (by a double) twice
// as many as a -- so each operation's result type carries that bound:
//
//   grow(Expansion<N>, double)            -> Expansion<N + 1>
//   sum(Expansion<M>, Expansion<N>)       -> Expansion<M + N>
//   fast_sum(Expansion<M>, Expansion<N>)  -> Expansion<M + N>
//   scale(Expansion<N>, double)           -> Expansion<2 * N>
//   product(Expansion<M>, Expansion<N>)   -> Expansion<2 * M * N>
//   negate(Expansion<N>)                  -> Expansion<N>
//
// and the operators +, -, * map onto these with the same capacities. The one
// way down is fit<K>(e), checked at runtime, for when a smaller bound is a
// theorem the types cannot see.
//
// A predicate written as a formula over these then has its sizes worked out
// by the compiler, and an expression that could overflow its storage does not
// compile: a sizing mistake is a type error, not a buffer overrun at runtime
// -- which is the bookkeeping Shewchuk's C code does by hand, in comments.
#pragma once

#include "kernel/eft.h"
#include <array>
#include <cassert>
#include <cstddef>

// An exact number, as at most N nonoverlapping doubles in increasing
// magnitude. See the top of this file for the invariant and the capacity.
template <std::size_t N> class Expansion
{
  public:
    static constexpr std::size_t capacity = N;

    Expansion() = default;       // zero: no components
    explicit Expansion(double a) // one component, or none if a == 0
    {
        append(a);
    }
    Expansion(TwoTerm t) // {t.lo, t.hi}, zeros dropped
    {
        append(t);
    }
    // Widening: the same value, in a larger capacity. For when two results of
    // different capacities have to be stored alike -- the terms of a sum
    // over a determinant's cofactors, say. Narrowing has no such constructor,
    // because it could not be done without checking at runtime.
    template <std::size_t M>
    Expansion(const Expansion<M> &e) // widening: requires M <= N
    {
        static_assert(
            M <= N,
            "an Expansion<M> only widens into an Expansion<N> with N >= M");
        for (std::size_t i = 0; i < e.size(); i++)
        {
            append(e[i]);
        }
    }

    std::size_t size() const
    {
        return count_;
    }

    // By value, on purpose: a component can be read but not written, so
    // nothing outside append() can break the invariant.
    double operator[](std::size_t i) const // smallest first
    {
        return c_[i];
    }
    // Adds `x` as the new largest component; does nothing if x is zero.
    //
    // For the operations in this file, which produce their components
    // smallest first and nonoverlapping, as the theorems guarantee. It
    // asserts what it can check cheaply -- that there is room, and that x is
    // larger than the last component -- and the tests keep the asserts on, so
    // an operation that breaks the invariant fails at the component that
    // breaks it.
    void append(double x)
    {
        if (x == 0)
        {
            return;
        }
        assert(count_ < N);
        assert(count_ == 0 || std::abs(x) > std::abs(c_[count_ - 1]));
        c_[count_] = x;
        count_++;
    }

    // A two-term from eft.h: its error first, then its rounded value.
    void append(TwoTerm term)
    {
        append(term.lo);
        append(term.hi);
    }

  private:
    std::array<double, N> c_{}; // components
    std::size_t count_ = 0;     // how many components
};

// e + b, exactly. Shewchuk's Grow-Expansion (Theorem 10), keeping no zeros.
//
// b is carried up through the components, smallest first. At each step
// two_sum(q, e[i]) splits q + e[i] into its rounded value, which becomes the
// new carry, and its error, which is final: everything still to come is
// larger than it and does not overlap it. So the errors come out in
// increasing order, and the carry left at the end is the largest component.
//
// two_sum, not fast_two_sum: the carry can be larger or smaller than e[i], so
// fast_two_sum's precondition cannot be promised.
template <std::size_t N> Expansion<N + 1> grow(const Expansion<N> &e, double b)
{
    Expansion<N + 1> res;
    double q = b;
    for (std::size_t i = 0; i < e.size(); i++)
    {
        TwoTerm term = two_sum(q, e[i]);
        q = term.hi;
        res.append(term.lo);
    }
    res.append(q);
    return res;
}

// a + b, exactly. Shewchuk's Expansion-Sum (Theorem 12), keeping no zeros.
//
// a is grown by each component of b in turn, smallest first: a sum of
// expansions is a sequence of grows. It cannot be written as a loop of grow()
// calls, though, because each grow returns a larger type and a loop cannot
// change a variable's type. So the grows run in place, in one buffer of the
// final capacity, each pass the same carry loop as grow().
//
// The buffer is a plain array rather than an Expansion because, during a
// pass, it does not satisfy the invariant: the new components are written
// back over the old ones from the start (k never passes i, so nothing is
// overwritten before it is read), and only once the pass ends is the prefix
// [0, n) an expansion again. The result is then appended into a real one.
//
// Each pass regrows from the bottom of the buffer. Shewchuk's version starts
// pass j at component j instead, since the components below it are already
// final; with zeros eliminated those positions are no longer fixed, so this
// simpler form is kept here. It is O(m * n) either way, which is what step 3,
// Fast-Expansion-Sum, improves on -- and this version is what the tests check
// that one against.
template <std::size_t M, std::size_t N>
Expansion<M + N> sum(const Expansion<M> &a, const Expansion<N> &b)
{
    std::array<double, N + M> buffer;
    for (std::size_t i = 0; i < a.size(); i++)
    {
        buffer[i] = a[i];
    }
    std::size_t n = a.size();
    for (std::size_t j = 0; j < b.size(); j++)
    {
        double q = b[j];
        std::size_t k = 0;
        for (std::size_t i = 0; i < n; i++)
        {
            TwoTerm term = two_sum(q, buffer[i]);
            q = term.hi;
            if (term.lo != 0)
            {
                buffer[k] = term.lo;
                k++;
            }
        }
        if (q != 0)
        {
            buffer[k] = q;
            k++;
        }
        n = k;
    }

    Expansion<M + N> res;
    for (std::size_t i = 0; i < n; i++)
    {
        res.append(buffer[i]);
    }
    return res;
}
// e * b, exactly, for a double b. Shewchuk's Scale-Expansion (Theorem 19),
// keeping no zeros.
//
// Each component times b is two_product's exact two parts, T (rounded) and t
// (its error), and those are folded into a running carry with two steps:
//
//   two_sum(carry, t)        t and the carry are of similar size, and either
//                            can be the larger: two_sum. Its error is final.
//   fast_two_sum(T, s)       T is the product of a larger component than any
//                            before it, so it dominates what has been carried
//                            -- Shewchuk's proof shows the precondition holds
//                            -- and the cheap version is safe. Its error is
//                            final too, and its rounded value is the new carry.
//
// Two outputs per component is where the capacity 2N comes from. The first
// component has no carry to fold into, so it is taken on its own before the
// loop. Scaling by zero is the empty expansion, returned at once rather than
// as a list of zero products.
//
// Exact within eft.h's range for products: the error terms must stay above
// the subnormals, so the components times b must not fall below 2^-969.
template <std::size_t N> Expansion<2 * N> scale(const Expansion<N> &e, double b)
{
    Expansion<2 * N> res;
    if (e.size() == 0 || b == 0)
    {
        return res;
    }
    TwoTerm term = two_product(e[0], b);
    res.append(term.lo);
    for (std::size_t i = 1; i < e.size(); i++)
    {
        TwoTerm t1 = two_product(e[i], b);
        TwoTerm t2 = two_sum(term.hi, t1.lo);
        res.append(t2.lo);
        term = fast_two_sum(t1.hi, t2.hi);
        res.append(term.lo);
    }
    res.append(term.hi);
    return res;
}

// The sign of the exact value: -1, 0 or +1. The sign of the largest
// component, which the rest cannot outweigh (see the top of this file); 0 for
// the empty expansion, the only representation zero has.
template <std::size_t N> int sign(const Expansion<N> &e)
{
    if (e.size() == 0)
    {
        return 0;
    }
    else if (e[e.size() - 1] > 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

// -e, exactly. Negating every component keeps the invariant: nonoverlapping
// is about which bits are set, not about signs.
template <std::size_t N> Expansion<N> negate(const Expansion<N> &e)
{
    Expansion<N> res;
    for (std::size_t i = 0; i < e.size(); i++)
    {
        res.append(-e[i]);
    }
    return res;
}

// The exact value, rounded to about a double: the components added smallest
// first. Shewchuk's "Estimate". Not exact, and not meant to be -- it is for
// reading an expansion, by a filter or a constructed position -- but close:
// the small components are added together first, while they still fit, so
// little is lost before the large ones arrive. The tests hold it within 2^-51
// of the exact value, relatively. (The largest component alone is no
// substitute; see the top of this file.)
template <std::size_t N> double approximate(const Expansion<N> &e)
{
    double res = 0.0;
    for (std::size_t i = 0; i < e.size(); i++)
    {
        res += e[i];
    }
    return res;
}

// a + b, exactly, in one pass. Shewchuk's Fast-Expansion-Sum (Theorem 13),
// keeping no zeros.
//
// sum() grows a by each component of b separately: one full carry pass per
// component, O(m * n). This merges the two expansions into one list ordered
// by magnitude -- the merge step of merge sort -- and runs a single carry
// pass up it, as grow() does: O(m + n). Every exact sum of two expansions in
// a predicate is this one.
//
// The merged list is not an expansion: a component of a and one of b can
// share bits. That one carry pass still produces a valid expansion is the
// content of the theorem, and it holds under two conditions:
//
//   round-to-nearest-even, as everything in eft.h assumes;
//
//   strongly nonoverlapping inputs (Shewchuk, section 2.3): nonoverlapping,
//   and where two neighbouring components touch -- one's highest bit just
//   below the other's lowest -- both are powers of two, and no component
//   touches both its neighbours. grow(), sum() and scale() produce more than
//   that (nonadjacent: no neighbours touch at all), so their results qualify.
//
// The result is strongly nonoverlapping too, but not necessarily nonadjacent:
// enough to be fed back into fast_sum, and to everything else here, which
// needs only nonoverlapping. tests/test_fast_sum.cpp checks exactly this
// contract, including on a tree of fast_sums whose every level eats the last.
//
// The carry pass starts with fast_two_sum: after the merge, g[1] is at least
// as large as g[0], so its precondition holds. From then on the carry gathers
// everything below and can outgrow the next component, so the order is
// unknown and it has to be two_sum. The merge buffer is a plain array, like
// sum()'s, because the merged list does not satisfy the invariant; only the
// result is built through append().
template <std::size_t M, std::size_t N>
Expansion<M + N> fast_sum(const Expansion<M> &a, const Expansion<N> &b)
{
    std::array<double, M + N> g{};
    std::size_t i = 0, j = 0, n = 0;
    // 1. merge: take whichever input's component is smaller in magnitude
    while (i < a.size() && j < b.size())
    {
        if (std::abs(a[i]) < std::abs(b[j]))
        {
            g[n] = a[i];
            i++;
        }
        else
        {
            g[n] = b[j];
            j++;
        }
        n++;
    }
    for (std::size_t ci = i; ci < a.size(); ci++)
    {
        g[n] = a[ci];
        n++;
    }
    for (std::size_t cj = j; cj < b.size(); cj++)
    {
        g[n] = b[cj];
        n++;
    }

    // 2. one carry pass, as in grow
    Expansion<M + N> res;
    if (n == 0)
    {
        return res;
    }
    double q = g[0];
    if (n >= 2)
    {
        TwoTerm term = fast_two_sum(g[1], q);
        q = term.hi;
        res.append(term.lo);
    }
    for (std::size_t k = 2; k < n && k < g.size(); k++)
    {
        TwoTerm term = two_sum(q, g[k]);
        q = term.hi;
        res.append(term.lo);
    }
    res.append(q);
    return res;
}

// e, moved into an Expansion<K> -- which may be smaller than e's capacity.
//
// The one narrowing, and checked at runtime: it asserts that e's components
// fit. Every capacity in this file is an upper bound worked out from the
// types, and the types only see the operations, not their meaning. Sometimes a
// tighter bound is a theorem -- in product() below, a running sum after j
// terms needs at most 2M(j + 1) components, whatever the types say -- and fit
// is how that is told to the compiler, with the assert standing in for the
// proof. The tests keep the assert on, so a wrong bound fails loudly.
//
// K is given, N is deduced: fit<6>(e). That is why K comes first: explicit
// template arguments fill the parameters from the left.
template <std::size_t K, std::size_t N> Expansion<K> fit(const Expansion<N> &e)
{
    Expansion<K> res;
    assert(e.size() <= K);
    for (std::size_t i = 0; i < e.size(); i++)
    {
        res.append(e[i]);
    }
    return res;
}

// a * b, exactly, for two expansions: the sum over b's components of
// scale(a, b[j]), accumulated with fast_sum.
//
// The accumulator would grow by a type with every term -- fast_sum returns a
// larger capacity than its inputs -- so each partial sum is brought back to
// the final capacity with fit. That is safe because the bound is real: each
// scale() gives at most 2M components and fast_sum never gives more than its
// inputs together, so after j + 1 terms there are at most 2M(j + 1) <= 2MN.
//
// fast_sum's precondition holds at every step: the accumulator is a fast_sum
// result (strongly nonoverlapping) and each scale() result is nonadjacent.
// Exact within scale()'s range: the components' products must stay above
// 2^-969 (see eft.h).
//
// Shewchuk's code has hand-written products of fixed sizes (Two_Two_Product
// and the like) that are faster for two-component factors. This one is
// general, and is what a predicate formula's * means; specialise it only if a
// profile says so.
template <std::size_t M, std::size_t N>
Expansion<2 * M * N> product(const Expansion<M> &a, const Expansion<N> &b)
{
    Expansion<2 * M * N> acc;
    for (std::size_t j = 0; j < b.size(); j++)
    {
        acc = fit<2 * M * N>(fast_sum(acc, scale(a, b[j])));
    }
    return acc;
}

// --- operators ---
//
// So that a predicate's formula can be written once, as mathematics, and
// evaluated with any number type: (ax - cx) * (by - cy) reads the same for
// doubles, where it rounds, and for expansions, where every operator below is
// one of the exact operations above and the result type -- its capacity --
// is worked out by the compiler (see orient2d.h). Later, the error-bounded
// double of the filter will be a third such type.
//
// Each operator is one exact operation, with that operation's capacity:
//
//   a + b   fast_sum(a, b)            Expansion<M + N>
//   a - b   fast_sum(a, negate(b))    Expansion<M + N>
//   -a      negate(a)                 Expansion<N>
//   a * b   product(a, b)             Expansion<2 * M * N>
//   a * d   scale(a, d), d a double   Expansion<2 * N>
//
// There is deliberately no mixing of an expansion and a double in + and -:
// a formula's inputs are made Expansion<1> first, so every operand is
// already an expansion, and the capacities stay predictable.
template <std::size_t N, std::size_t M>
Expansion<M + N> operator+(const Expansion<N> &a, const Expansion<M> &b)
{
    return fast_sum(a, b);
}

template <std::size_t N, std::size_t M>
Expansion<M + N> operator-(const Expansion<N> &a, const Expansion<M> &b)
{
    return fast_sum(a, negate(b));
}

template <std::size_t N, std::size_t M>
Expansion<2 * M * N> operator*(const Expansion<N> &a, const Expansion<M> &b)
{
    return product(a, b);
}

template <std::size_t N>
Expansion<2 * N> operator*(const Expansion<N> &a, double d)
{
    return scale(a, d);
}

template <std::size_t N> Expansion<N> operator-(const Expansion<N> &a)
{
    return negate(a);
}
