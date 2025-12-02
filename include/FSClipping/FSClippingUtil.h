#ifndef FSCLIPPINGUTIL_HPP
#define FSCLIPPINGUTIL_HPP

#include "FSClippingFace.h"

_FS_BEGIN_NAMESPACE

namespace FSClippingUtil {

/* -----------------------------------------------------------------
   AreFacesIdentical
   ----------------------------------------------------------------- */
inline bool AreFacesIdentical(const FSClippingFace &f1,
                              const FSClippingFace &f2,
                              FS_floatT tol = static_cast<FS_floatT>(1e-12)) {

  // improve this algorithm to have nlog(n) cost

  const auto &v1 = f1.vertices();
  const auto &v2 = f2.vertices();

  if (v1.size() != v2.size())
    return false;

  const FS_intT n = static_cast<FS_intT>(v1.size());
  std::vector<bool> used(n, false);

  for (FS_intT i = 0; i < n; ++i) {
    const auto &p = v1[i];
    bool matched = false;
    for (FS_intT j = 0; j < n; ++j) {
      if (used[j])
        continue;
      if ((p - v2[j]).L2Norm() < tol) {
        used[j] = true;
        matched = true;
        break;
      }
    }
    if (!matched)
      return false;
  }
  return true;
}

/* -----------------------------------------------------------------
   IsPointInsideConvexPolygon
   ----------------------------------------------------------------- */
inline bool
IsPointInsideConvexPolygon(const FSVec2 &P, const std::vector<FSVec2> &poly,
                           FS_floatT tol = static_cast<FS_floatT>(1e-12)) {
  const FS_intT n = static_cast<FS_intT>(poly.size());
  for (FS_intT i = 0; i < n; ++i) {
    const FSVec2 &A = poly[i];
    const FSVec2 &B = poly[(i + 1) % n];
    FSVec2 edge = B - A;
    // normale intérieure (orientation CCW)
    FSVec2 normal(-edge[1], edge[0]);
    FS_floatT d = (P - A).InnerProduct(normal);
    if (d < -tol)
      return false; // à l’extérieur
  }
  return true;
}

/* -----------------------------------------------------------------
   IsFaceInsideTheOther
   ----------------------------------------------------------------- */
inline bool
IsFaceInsideTheOther(const std::vector<FSVec2> &f1,
                     const std::vector<FSVec2> &f2,
                     FS_floatT tol = static_cast<FS_floatT>(1e-12)) {
  for (const FSVec2 &P : f1) {
    if (!IsPointInsideConvexPolygon(P, f2, tol))
      return false;
  }
  return true;
}

} // namespace FSClippingUtil

_FS_END_NAMESPACE
#endif // FSCLIPPINGUTIL_HPP