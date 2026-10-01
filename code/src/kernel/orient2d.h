// orient2d: which side of a line a point is on -- the first predicate, and
// the pattern every later one follows.
//
// orient2d(a, b, c) is the sign of the determinant
//
//   | ax - cx   ay - cy |
//   | bx - cx   by - cy |  =  (ax - cx)(by - cy) - (ay - cy)(bx - cx),
//
// twice the signed area of the triangle abc: positive when a, b, c turn
// counterclockwise, negative when they turn clockwise, zero when they are
// collinear. Shewchuk's convention.
//
// --- one formula, several number types ---
//
// The determinant is written once, in orient2d_det, as a template. What it
// computes depends on the type it is given:
//
//   double         the plain formula: fast, and rounded. Right almost
//                  everywhere, and wrong where it matters -- a few ulps from
//                  the line, where the rounding error is larger than the
//                  determinant. tests/test_orient2d.cpp shows it failing on
//                  the grid of Kettner et al.'s "Classroom examples".
//
//   Expansion<1>   the same text, exact: each - and * is an exact operation
//                  of expansion.h, and the compiler works out the sizes --
//                  differences Expansion<2>, products Expansion<8>, the
//                  determinant Expansion<16>. No component is ever rounded,
//                  so the sign is the true sign, and collinear points give
//                  exactly zero.
//
// Step 5 adds a third type: a double that carries a bound on its own rounding
// error -- the filter. orient2d will try that first and fall back to the
// expansion only when the bound cannot decide, which on ordinary input is
// almost never. The formula will not change; that is the point of writing it
// once.
//
// --- why the differences are taken against c ---
//
// Translating the three points so that c is the origin is Shewchuk's form of
// the determinant. The untranslated one, the 3 x 3 determinant with a column
// of ones, expands into six products of coordinates; exactly, the two forms
// cost about the same. The difference is in the double evaluation, which is
// what the filter will rely on: there, each difference is computed first and
// rounds by at most half an ulp of itself, so the rounding error of the whole
// determinant is proportional to the *differences* -- the size of the
// triangle -- and not to the coordinates. Three points near (1e6, 1e6) that
// are a millimetre apart then give an error bound on the scale of their own
// tiny area, and the filter can still decide them; the untranslated form's
// bound would be on the scale of 1e6 * 1e6, and decide almost nothing there.
//
// --- range ---
//
// Exact as long as no product overflows or underflows (kernel/eft.h): every
// coordinate zero, or of magnitude between about 2^-400 and 2^500, is
// comfortably inside. Mesh and seed coordinates are, and the domain and seed
// set will check it where they are built.
#pragma once
#include "kernel/expansion.h"
#include <type_traits>

// The determinant, for any number type with -, * and copy. The return type is
// whatever the formula produces: a double for doubles, an Expansion<16> for
// Expansion<1> inputs.
template <class T>
auto orient2d_det(const T &ax, const T &ay, const T &bx, const T &by,
                  const T &cx, const T &cy)
{
    return (ax - cx) * (by - cy) - (ay - cy) * (bx - cx);
}

// The exact sign: +1 counterclockwise, -1 clockwise, 0 collinear.
//
// Each coordinate is made an Expansion<1> -- a double is exactly its own
// value, so nothing is lost -- and the same template evaluates the
// determinant exactly. Always the exact path, for now; step 5 puts the filter
// in front of it. inline because it is a plain function defined in a header:
// every file that includes this one gets a copy, and inline tells the linker
// they are all the same function.
inline int orient2d(double ax, double ay, double bx, double by, double cx,
                    double cy)
{
    using E = Expansion<1>;
    const auto det = orient2d_det(E(ax), E(ay), E(bx), E(by), E(cx), E(cy));
    static_assert(std::is_same_v<decltype(det), const Expansion<16>>);
    return sign(det);
}
