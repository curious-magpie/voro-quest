// Error-free transformations: one rounded floating-point operation, turned into
// two doubles whose sum is the exact result.
//
// Every double operation rounds. a + b computed in double is fl(a + b), the
// true sum rounded to the nearest double, and the difference is lost. An
// error-free transformation (EFT) keeps it:
//
//   two_sum(a, b)      = {hi, lo}  with  hi = fl(a + b)  and  hi + lo = a + b
//   two_product(a, b)  = {hi, lo}  with  hi = fl(a * b)  and  hi + lo = a * b
//
// where "=" on the right is exact, in the real numbers. Nothing is
// approximated: the rounding error of the operation is itself a double, and it
// is computed exactly, from doubles, with nothing but more rounded operations.
// Everything exact in RVD -- the expansions of expansion.h, and through them
// every predicate -- is built from these five functions.
//
// The name is from Ogita, Rump and Oishi (2005). The algorithms are older:
// Dekker (1971) and Knuth (TAOCP, vol. 2). The proofs this file leans on are in
//
//   J. R. Shewchuk, Adaptive precision floating-point arithmetic and fast
//   robust geometric predicates, Discrete & Computational Geometry 18, 1997,
//   section 2,
//
// whose theorem numbers are quoted below.
//
// --- what "exact" rests on ---
//
// The guarantees hold only under the arithmetic they were proved for:
//
//   IEEE 754 binary64 with round-to-nearest-even, the default everywhere.
//
//   No compiler rewriting. Each function is a sequence of operations whose
//   order is the algorithm; in real arithmetic most of them cancel out (lo is
//   always 0 in the reals), so an optimiser allowed to use real-number algebra
//   would delete them. Never build this with -ffast-math. And -ffp-contract=off
//   (set for the whole project in CMakeLists.txt) stops the compiler fusing a
//   multiply and an add into one fma, which rounds once instead of twice and
//   breaks the proofs. An explicit std::fma, as in two_product, is still
//   honoured: it is asked for, not invented.
//
//   No extended precision. On the old x87 FPU, intermediates were kept in 80
//   bits and rounded again on store ("double rounding"), which breaks these
//   algorithms. x86-64 and ARM compute doubles in SSE2/NEON registers, so this
//   no longer arises, but it is why old robust-predicate code is full of
//   volatile.
//
//   Finite inputs, and these ranges:
//     two_sum, fast_two_sum   exact unless the sum overflows. Gradual
//                             underflow keeps sums exact even among subnormals.
//     two_product             exact while the error term stays out of the
//                             subnormal range: |a * b| >= 2^-969 (that is,
//                             2^(-1022 + 53)), or a * b == 0.
//     split                   overflows when |a| > 2^996, about 6.7e299:
//                             a * (2^27 + 1) must stay finite.
//     two_product_dekker      both of the above, for a and for b.
//   Mesh and seed coordinates are nowhere near these limits, but a predicate
//   of degree 6 multiplies six of them, so the limits are checked where a
//   domain or a seed set is built (design, section 6), not here.
//
// --- the facts the proofs use ---
//
// Writing ulp(x) for the gap between x and the next double away from zero:
//
//   (F1) The rounding error of a sum is a double. For doubles a and b, the
//        exact a + b - fl(a + b) is itself representable, and at most half an
//        ulp of fl(a + b). (Round-to-nearest cannot lose more than half an ulp,
//        and the lost part is made of bits a and b already had.)
//
//   (F2) Sterbenz's lemma. If b/2 <= a <= 2b, then a - b is computed exactly:
//        the result needs no more bits than the operands, so nothing rounds.
//
//   (F3) The rounding error of a product is a double, provided it does not
//        underflow: a * b - fl(a * b) has at most 53 significant bits, because
//        the exact product of two 53-bit significands has at most 106.
//
// --- "nonoverlapping", which expansions depend on ---
//
// Two doubles x and y are nonoverlapping (Shewchuk, section 2.2) when the
// lowest set bit of the larger is above the highest set bit of the smaller:
// they cover disjoint ranges of binary digits, as 1000 and 11 do. Every {hi,
// lo} returned here is nonoverlapping, with the stronger property that
// |lo| <= ulp(hi) / 2, so hi + lo rounds back to hi. That is what lets an
// expansion -- a list of such parts -- be read off: its sign is the sign of its
// largest part, and its largest part is a good approximation of the whole.
// tests/test_eft.cpp checks exactly this, for every function here.
#pragma once

constexpr double kSplitter = 134217729.0; // 2^27 + 1

#include <cassert>
#include <cmath>

// An exact value as two parts: hi, the rounded result, and lo, what the
// rounding lost. hi + lo is the exact result; hi and lo do not overlap.
struct TwoTerm
{
    double hi;
    double lo;
};

// Shewchuk's Theorem 6 (Dekker): if |a| >= |b|, then {hi, lo} is a
// nonoverlapping expansion with hi + lo = a + b exactly.
//
// Three operations, the cheapest EFT there is, and the precondition is the
// price. Why it works: hi = fl(a + b) is a + b rounded, and because b is the
// smaller term, hi - a is computed exactly (it is close to b, and it needs no
// more bits than hi and a have in common -- the argument is a case analysis
// on the signs, of the kind F2 settles). So bv is *exactly* the part of b that
// made it into hi, and b - bv is the part that did not: the rounding error,
// which is a double by F1, so that subtraction is exact too.
//
// Without the precondition, hi - a can round, bv is no longer the part of b
// that made it in, and lo is simply wrong. Use two_sum when the order is not
// known; expansion code uses this where its own invariants guarantee it.
inline TwoTerm fast_two_sum(double a, double b)
{
    assert(std::abs(a) >= std::abs(b));
    TwoTerm res;
    res.hi = a + b;
    double bv = res.hi - a; // exact: the part of b that made it into hi
    res.lo = b - bv;        // exact: the part of b that did not
    return res;
}

// Shewchuk's Theorem 7 (Knuth): for any a and b, {hi, lo} is a nonoverlapping
// expansion with hi + lo = a + b exactly.
//
// fast_two_sum without the precondition, for three more operations. Neither
// operand is assumed to be the larger, so both are "recovered" from hi:
//
//   bv = hi - a     the part of b that hi carries ("b virtual")
//   av = hi - bv    the part of a that hi carries ("a virtual")
//   br = b - bv     what was lost of b ("b roundoff")
//   ar = a - av     what was lost of a ("a roundoff")
//   lo = ar + br    the whole rounding error
//
// Individually, bv and av may themselves be rounded when |b| > |a|. The
// theorem's content is that the errors made in recovering them cancel
// exactly in ar + br, so lo is the exact error all the same. Note av uses bv,
// not b: subtracting b itself is not exact in general, and the cancellation
// argument needs the recovered value.
inline TwoTerm two_sum(double a, double b)
{
    TwoTerm res;
    res.hi = a + b;
    double bv = res.hi - a;
    double av = res.hi - bv;
    double ar = a - av;
    double br = b - bv;
    res.lo = ar + br;
    return res;
}

// For any a and b in range (see the top of this file), {hi, lo} is a
// nonoverlapping expansion with hi + lo = a * b exactly.
//
// The modern way, and two operations instead of Dekker's seventeen. std::fma
// computes a * b - hi with a single rounding. That exact value is a double by
// F3, and rounding a double changes nothing, so the fma returns it exactly.
// hi is computed once and reused, so hi and lo describe the same rounded
// product.
//
// std::fma is exact by the C++ standard on every platform, but only fast where
// the CPU has an fma instruction (every x86-64 since Haswell, every 64-bit
// ARM). Elsewhere it is emulated in software, slowly -- which is why
// two_product_dekker is kept, and tested to agree bit for bit.
inline TwoTerm two_product(double a, double b)
{
    TwoTerm res;
    res.hi = a * b;
    res.lo = std::fma(a, b, -res.hi);
    return res;
}

// Shewchuk's Theorem 17 (Dekker, after Veltkamp): splits a into hi + lo = a
// exactly, each half with at most 26 significant bits.
//
// The point is two_product_dekker: two 26-bit numbers multiply to at most 52
// bits, which fits in a double's 53, so every product of halves is exact.
//
// How it works. Multiplying by 2^27 + 1 gives c = a * 2^27 + a, rounded. In
// c - a, the copy of a mostly cancels, and what is left is a * 2^27 with its
// low bits rounded away; subtracting that from c leaves hi: a rounded to its
// top 53 - 27 = 26 bits. Then lo = a - hi is exact (F2 territory: hi is
// within half its last bit of a). Shewchuk's theorem, for p-bit doubles and
// splitter 2^s + 1, gives hi p - s bits and lo s - 1 bits; with p = 53 and
// s = 27 both have 26. lo holds 26 bits where 27 seem needed because its sign
// carries one: rounding to nearest leaves a remainder that may be negative.
//
// The constant is (1 << 27) + 1 = 134217729 = 2^27 + 1, with 27 = ceil(53/2):
// the smallest split point that leaves both halves at most 26 bits.
inline TwoTerm split(double a)
{
    TwoTerm res;
    double c = a * kSplitter;
    double abig = c - a;
    res.hi = c - abig;
    res.lo = a - res.hi;
    return res;
}

// Shewchuk's Theorem 18 (Dekker): for any a and b in range (see the top of
// this file), {hi, lo} is a nonoverlapping expansion with hi + lo = a * b
// exactly -- the same result as two_product, without fma.
//
// Split both operands, so that a * b is the exact sum of four products of
// halves, each exact in a double:
//
//   a * b = ahi * bhi + alo * bhi + ahi * blo + alo * blo
//
// and the error a * b - hi is what remains when hi is taken away. hi is close
// to ahi * bhi, the largest of the four, so hi - ahi * bhi is exact (F2), and
// each further subtraction, largest product first, peels off one term without
// rounding: the remainder keeps getting smaller and keeps fitting. The last
// line flips the sign, so that lo = a * b - hi rather than its negative.
//
// The order matters: subtract the products in any other order and the
// intermediate remainders can need more than 53 bits.
inline TwoTerm two_product_dekker(double a, double b)
{
    TwoTerm res;
    res.hi = a * b;
    TwoTerm asplit = split(a);
    TwoTerm bsplit = split(b);
    double err1 = res.hi - (asplit.hi * bsplit.hi);
    double err2 = err1 - (asplit.lo * bsplit.hi);
    double err3 = err2 - (asplit.hi * bsplit.lo);
    res.lo = (asplit.lo * bsplit.lo) - err3;
    return res;
}
