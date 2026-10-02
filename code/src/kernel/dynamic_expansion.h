// Expansions whose length is decided at runtime.
//
// The same exact numbers as expansion.h -- a list of nonzero, nonoverlapping
// doubles in increasing magnitude whose sum is the value -- and the same
// algorithms, Shewchuk's. What changes is where the capacity lives.
// Expansion<N> carries its bound in the type, which is right for orient2d and
// orient3d: the compiler works the sizes out, and nothing is ever allocated.
// But the bound multiplies with every product, and the clipping predicates of
// step 6c are polynomials of degree 4 to 6: side2 would need an
// Expansion<10000> and side3 about three million components, far more than
// the stack holds. Those are worst cases of the formula; the real length of
// an expansion is set by how many bits its value has, which is usually far
// less, and known only once the value is computed. So here the length is a
// runtime number.
//
// --- the storage ---
//
// A DynamicExpansion keeps its first 32 components in an array inside the
// object, and moves all of them to a std::vector once a 33rd arrives. The
// expectation is that most exact evaluations are short -- the filter of
// bounded.h only falls back to the exact path near a degeneracy -- and those
// never touch the heap. The long ones pay for one allocation, which is cheap
// next to the arithmetic that made them long.
//
// All the components are in one place or the other, never split, and nothing
// points into the inline array. That keeps operator[] to one branch, and it
// is what makes the copy and move the compiler writes correct: a pointer into
// inline_ would, in a copy, still point into the original's array. What the
// default move does not do is leave the source readable: after
// `D b = std::move(a)` on a heap-stored a, a keeps its count while its vector
// is (in practice) empty. As with any moved-from object, assign to it before
// reading it.
//
// --- the operations ---
//
// Free functions with the names of expansion.h's, so a predicate's formula
// reads the same with either type. Two differ, because the capacity is no
// longer in the type: product needs no fit, and there is compress, which
// Expansion<N> does not have. Each is `inline` because it is defined in a
// header: a non-template function defined in a header that two .cpp files
// include would otherwise be defined twice, and the link would fail.
#pragma once
#include <array>
#include <cassert>
#include <vector>
#include <cmath>
#include <cstddef>
#include "kernel/eft.h"
#include "kernel/expansion.h"

// An exact number, as nonoverlapping doubles in increasing magnitude, as many
// as it needs. See the top of this file for the storage, and the top of
// expansion.h for the invariant, which is the same.
class DynamicExpansion
{
  public:
    // How many components fit without allocating. A first guess, to be
    // revisited with a profile of the side predicates.
    static constexpr std::size_t kInline = 32;

    DynamicExpansion() = default;       // zero: no components
    explicit DynamicExpansion(double a) // one component, or none if a == 0
    {
        append(a);
    }
    DynamicExpansion(TwoTerm t) // {t.lo, t.hi}, zeros dropped
    {
        append(t.lo);
        append(t.hi);
    }
    // From a fixed-capacity expansion: the same components, so that a value
    // computed with Expansion<N> can continue here.
    template <std::size_t N>
    explicit DynamicExpansion(const Expansion<N> &e) // append each e[i]
    {
        for (std::size_t i = 0; i < e.size(); i++)
        {
            append(e[i]);
        }
    }

    std::size_t size() const
    {
        return count_;
    }
    // By value, as in Expansion: nothing outside append() can write a
    // component, so nothing else can break the invariant. Where the component
    // is read from depends only on the count: inline up to kInline, the heap
    // above, matching where append() put them all.
    double operator[](std::size_t i) const // smallest first
    {
        assert(i < count_);
        if (count_ <= kInline)
        {
            return inline_[i];
        }
        else
        {
            return heap_[i];
        }
    }
    // Adds x as the new largest component; does nothing if x is zero. The
    // same checks as Expansion::append, minus the capacity, which is
    // unbounded.
    //
    // The 33rd component is the one that moves storage: the 32 inline ones
    // are copied into the vector first, then x is added after them, so that
    // from then on every component is on the heap. The inline array is left
    // as it was and is no longer read.
    void append(double x)
    {
        if (x == 0)
        {
            return;
        }
        assert(count_ == 0 || std::abs(x) > std::abs((*this)[count_ - 1]));
        if (count_ < kInline)
        {
            inline_[count_] = x;
        }
        else if (count_ == kInline)
        {
            heap_.assign(std::begin(inline_), std::end(inline_));
            heap_.push_back(x);
        }
        else
        {
            heap_.push_back(x);
        }
        count_++;
    }

  private:
    std::array<double, kInline> inline_{}; // components while count_ <= 32
    std::vector<double> heap_;             // all components once count_ > 32
    std::size_t count_ = 0;                // how many components
};

// e + b, exactly. Grow-Expansion, as in expansion.h: b is carried up through
// the components, and each two_sum's error is final.
inline DynamicExpansion grow(const DynamicExpansion &e, double b)
{
    DynamicExpansion res;
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

// a + b, exactly. Fast-Expansion-Sum, as in expansion.h, with the same
// contract: strongly nonoverlapping inputs give a nonoverlapping result. Not
// a strongly nonoverlapping one, whatever Theorem 13 says (see expansion.h),
// so the operators sum with linear_sum below instead.
//
// The difference is that the merged list is never stored. The carry pass
// reads it once, in order, so next() hands out the merge one component at a
// time: whichever of a[i] and b[j] is smaller in magnitude, or what is left of
// one input once the other is used up. A buffer of the merged list would be a
// heap allocation on every call, and a sum runs once per component of b in
// every product.
//
// next() is a lambda: a small function defined inside this one. [&] lets it
// use the surrounding variables by reference, so its i++ and j++ advance the
// real indices. It is called exactly n times, once per component, so neither
// input is read past its end.
inline DynamicExpansion fast_sum(const DynamicExpansion &a,
                                 const DynamicExpansion &b)
{
    const std::size_t n = a.size() + b.size();
    std::size_t i = 0, j = 0;
    auto next = [&]()
    {
        if (j == b.size() || (i < a.size() && std::abs(a[i]) < std::abs(b[j])))
        {
            return a[i++];
        }
        return b[j++];
    };

    DynamicExpansion res;
    if (n == 0)
    {
        return res;
    }
    // The first two merged components come out in order of magnitude, so
    // fast_two_sum's |a| >= |b| holds for them. After that the carry can
    // outgrow the next component, so it has to be two_sum.
    double q = next();
    if (n >= 2)
    {
        TwoTerm term = fast_two_sum(next(), q);
        q = term.hi;
        res.append(term.lo);
    }
    for (std::size_t k = 2; k < n; k++)
    {
        TwoTerm term = two_sum(q, next());
        q = term.hi;
        res.append(term.lo);
    }
    res.append(q);
    return res;
}

// a + b, exactly. Linear-Expansion-Sum (Theorem 24), as in expansion.h,
// with the same contract -- nonoverlapping inputs, which every
// DynamicExpansion is, give a nonoverlapping result -- and the same
// components: the merge breaks ties the same way, and the tests compare the
// two versions component for component. The sum the operators and product()
// use.
//
// The merge is fast_sum's next(), so nothing is buffered. It must be called
// exactly n times, and never twice in one expression: in
// fast_two_sum(next(), next()) C++ does not say which argument is evaluated
// first, so which component lands in which parameter would be unknown. So
// g[0] is read into g0 first, and each later component into gk. With one
// component in all, g0 is the sum.
inline DynamicExpansion linear_sum(const DynamicExpansion &a,
                                   const DynamicExpansion &b)
{
    const std::size_t n = a.size() + b.size();
    std::size_t i = 0, j = 0;
    auto next = [&]()
    {
        if (j == b.size() || (i < a.size() && std::abs(a[i]) < std::abs(b[j])))
        {
            return a[i++];
        }
        return b[j++];
    };

    DynamicExpansion res;
    if (n == 0)
    {
        return res;
    }

    double g0 = next();
    if (n == 1)
    {
        res.append(g0);
        return res;
    }

    TwoTerm t0 = fast_two_sum(next(), g0); // t0.hi is Q, t0.lo is q
    for (std::size_t k = 2; k < n; k++)
    {
        double gk = next();
        TwoTerm t1 = fast_two_sum(gk, t0.lo); // g[k] takes in q
        res.append(t1.lo);                    // final: below all that follows
        TwoTerm t2 = two_sum(t0.hi, t1.hi);   // the rest goes into Q
        t0 = t2;
    }
    res.append(t0.lo); // q, then Q: smallest first
    res.append(t0.hi);
    return res;
}

// e * b, exactly, for a double b. Scale-Expansion, as in expansion.h: each
// component's two_product is folded into the carry with a two_sum and a
// fast_two_sum, two final components per step. Scaling by zero is zero.
//
// Exact within eft.h's range for products: the components times b must not
// fall below 2^-969.
inline DynamicExpansion scale(const DynamicExpansion &e, double b)
{
    DynamicExpansion res;
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

// a * b, exactly, for two expansions: the sum over b's components of
// scale(a, b[j]), accumulated with linear_sum, as in expansion.h, so every
// step is a theorem (19, then 24). Without the fit: the accumulator is one
// type whatever its length, so the loop needs no help from the types.
//
// Not compressed; a * b is, see the operators below. The partial products
// are added one after another, so the accumulator is long for most of the
// loop; adding them in a balanced tree (pairs, then pairs of pairs) would do
// less work for long factors. Left as is until a profile of the side
// predicates says the exact path matters.
inline DynamicExpansion product(const DynamicExpansion &a,
                                const DynamicExpansion &b)
{
    DynamicExpansion acc;
    for (std::size_t j = 0; j < b.size(); j++)
    {
        acc = linear_sum(acc, scale(a, b[j]));
    }
    return acc;
}

// -e, exactly. Signs do not affect which bits are set, so the invariant holds.
inline DynamicExpansion negate(const DynamicExpansion &e)
{
    DynamicExpansion res;
    for (std::size_t i = 0; i < e.size(); i++)
    {
        res.append(-e[i]);
    }
    return res;
}

// The sign of the exact value: -1, 0 or +1. The sign of the largest
// component, which the rest cannot outweigh; 0 for the empty expansion.
inline int sign(const DynamicExpansion &e)
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

// The exact value, rounded to about a double: the components added smallest
// first, as in expansion.h. Close, not exact.
inline double approximate(const DynamicExpansion &e)
{
    double res = 0.0;
    for (std::size_t i = 0; i < e.size(); i++)
    {
        res += e[i];
    }
    return res;
}

// e, exactly, in fewer components. Shewchuk's Compress (section 2.8,
// Theorem 23).
//
// Products are where expansions get long: a product of an m- and an
// n-component expansion can have up to 2mn components, many carrying only a
// few bits that a neighbour could absorb. Compress merges neighbours whose
// sum is exactly a double, and gives two things back:
//
//   a nonadjacent result -- more than any operation here needs, fast_sum
//   included -- never longer than e;
//
//   a largest component within an ulp of the whole value. An expansion in
//   general does not have this ([1023.9, 1024] is worth twice its largest
//   component); a compressed one does, which makes it cheap to read.
//
// It does not promise the shortest expansion, and compressing again can
// change the length. Don't rely on either.
//
// Two passes over e:
//
//   pass 1, largest to smallest: Q gathers the components from the top down.
//   As long as fast_two_sum(Q, e[k]) has no error, e[k] is absorbed into Q.
//   When there is an error, Q is a finished component and is stored, from
//   the top of g downwards, and the error carries on as the new Q. At the
//   end, g[bottom .. m-1] holds the gathered components, largest at the top.
//
//   pass 2, smallest to largest: a second chance to merge, from below. Each
//   g[k] takes in the carry from beneath it, as in grow(); only the errors
//   that remain become components, and a zero error is dropped by append().
//
// Both passes use fast_two_sum, not two_sum: in pass 1 the running Q is at
// least as large as the next component below it, and in pass 2 g[k] is at
// least as large as the carry. Shewchuk's proof shows both, and the assert
// in fast_two_sum checks them in the tests.
//
// g holds up to m doubles and is written from the top, so it cannot be a
// DynamicExpansion, which only grows from the bottom. It is a plain buffer:
// an array on the stack for e of up to 32 components, a vector above. Both
// are declared at function scope, so that they live until the end -- a
// variable declared inside an if lives only until that if's closing brace --
// and g points at whichever is in use. An empty vector does not allocate, so
// short inputs never touch the heap.
//
// The pass 1 loop runs k from m-2 down to 0 with an unsigned k. The usual
// `k >= 0` is always true for an unsigned type, and k-- past 0 wraps around
// to a huge index; testing `k-- > 0` checks before decrementing, so the body
// sees m-2, ..., 0 and the loop stops after 0. For m == 1 it runs zero times.
inline DynamicExpansion compress(const DynamicExpansion &e)
{
    DynamicExpansion res;
    if (e.size() == 0)
    {
        return res;
    }
    std::size_t m = e.size();

    // pass 1 (downward, largest to smallest)
    std::array<double, DynamicExpansion::kInline> small;
    std::vector<double> big; // empty: this allocates nothing
    if (m > DynamicExpansion::kInline)
    {
        big.resize(m); // only long inputs pay for an allocation
    }
    double *g = (m <= DynamicExpansion::kInline) ? small.data() : big.data();
    std::size_t bottom = m - 1;
    double Q = e[m - 1];
    for (std::size_t k = m - 1; k-- > 0;)
    {
        TwoTerm term = fast_two_sum(Q, e[k]);
        Q = term.hi;
        if (term.lo != 0)
        {
            g[bottom] = Q;
            bottom--;
            Q = term.lo;
        }
    }
    g[bottom] = Q;

    // pass 2 (upward over g[bottom ... m-1])
    for (std::size_t k = bottom + 1; k < m; k++)
    {
        TwoTerm term = fast_two_sum(g[k], Q);
        Q = term.hi;
        res.append(term.lo);
    }
    res.append(Q);
    return res;
}

// --- operators ---
//
// As in expansion.h, so that a predicate's formula is written once and
// evaluated with doubles, the filter's bounded type, or exactly here:
//
//   a + b   linear_sum(a, b)
//   a - b   linear_sum(a, negate(b))
//   -a      negate(a)
//   a * b   compress(product(a, b))
//   a * d   scale(a, d), d a double
//
// The one difference is that * between expansions compresses. A product
// costs at least the product of its factors' lengths, so in a degree-6
// formula uncompressed factors would make every level slower than the one
// before; compressing costs one pass over the result and keeps each level
// short. Sums are not compressed: a sum is at most as long as its inputs
// together, and it is products that multiply lengths. Nor are the factors:
// with linear_sum inside product(), any nonoverlapping factor is a valid
// input, so compressing one would only be for speed -- and the one that would
// matter is b, whose length the cost grows with roughly quadratically. Left
// to a profile.
inline DynamicExpansion operator+(const DynamicExpansion &a,
                                  const DynamicExpansion &b)
{
    return linear_sum(a, b);
}

inline DynamicExpansion operator-(const DynamicExpansion &a,
                                  const DynamicExpansion &b)
{
    return linear_sum(a, negate(b));
}

inline DynamicExpansion operator*(const DynamicExpansion &a,
                                  const DynamicExpansion &b)
{
    return compress(product(a, b));
}

inline DynamicExpansion operator*(const DynamicExpansion &a, double d)
{
    return scale(a, d);
}

inline DynamicExpansion operator-(const DynamicExpansion &a)
{
    return negate(a);
}
