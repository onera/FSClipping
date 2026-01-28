#ifndef FSCLIPPINGUTIL_H
#define FSCLIPPINGUTIL_H

#include "FSPlaneSurface.h"

#include "FSClipping/FSClippingFace.h"

_FS_BEGIN_NAMESPACE

namespace FSClippingUtil {

/* -----------------------------------------------------------------
   PolygonSignedArea
   ----------------------------------------------------------------- */
inline FS_floatT PolygonSignedArea(const std::vector<FSVec2>& poly)
{
  const FS_intT n = static_cast<FS_intT>(poly.size());
  if(n < 3)
    return static_cast<FS_floatT>(0);

  FS_floatT area = 0;

  for(FS_intT i = 0; i < n; ++i) {
    const FSVec2& p = poly[i];
    const FSVec2& q = poly[(i + 1) % n];
    area += p[0] * q[1] - q[0] * p[1];
  }

  return static_cast<FS_floatT>(0.5) * area;
}

/* -----------------------------------------------------------------
   EnsureSameOrientation
   ----------------------------------------------------------------- */
inline void EnsureSameOrientation(std::vector<FSVec2>& poly,
                                  FS_floatT refSignedArea)
{
  FS_floatT area = PolygonSignedArea(poly);

  if(area * refSignedArea < 0.0) {
    std::reverse(poly.begin(), poly.end());
  }
}

/* -----------------------------------------------------------------
   IsPointInsideEdge
   ----------------------------------------------------------------- */
inline bool IsPointInsideEdge(const FSVec2& P,
                              const FSVec2& A,
                              const FSVec2& B,
                              FS_floatT tol = static_cast<FS_floatT>(1e-12))
{
  // (B-A) × (p-A) >= 0  →  p is on the left of AB
  FS_floatT cross =
    (B[0] - A[0]) * (P[1] - A[1]) - (B[1] - A[1]) * (P[0] - A[0]);
  return cross >= -tol;
}

/* -----------------------------------------------------------------
   ComputeIntersection
   ----------------------------------------------------------------- */
inline FSVec2 ComputeIntersection(const FSVec2& P1,
                                  const FSVec2& P2,
                                  const FSVec2& A,
                                  const FSVec2& B)
{
  FSVec2 U = P2 - P1; // direction of subject segment
  FSVec2 V = B - A;   // direction of clipping edge

  FS_floatT denom = U[0] * V[1] - U[1] * V[0];
  if(std::abs(denom) < static_cast<FS_floatT>(1e-12))
    return P1; // Parallel – return one endpoint as fallback

  FSVec2 AP = A - P1;
  FS_floatT t = (AP[0] * V[1] - AP[1] * V[0]) / denom;
  return P1 + U * t;
}

/* -----------------------------------------------------------------
   AreFacesIdentical
   ----------------------------------------------------------------- */
bool AreFacesIdentical(const FSClippingFace& f1,
                       const FSClippingFace& f2,
                       FS_floatT tol = static_cast<FS_floatT>(1e-12));

/* -----------------------------------------------------------------
   ClipConvexPolygon - Sutherland-Hodgman Algorithm
   For illustration : https://github.com/mhdadk/sutherland-hodgman
   ----------------------------------------------------------------- */
bool ClipConvexPolygon(const std::vector<FSVec2>& subject,
                       const std::vector<FSVec2>& clipper,
                       std::vector<FSVec2>& output,
                       FS_floatT tol = static_cast<FS_floatT>(1e-12));

/* -----------------------------------------------------------------
   ProjectFaceInPlane
   ----------------------------------------------------------------- */
std::vector<FSVec2> ProjectFaceInPlane(const FSClippingFace& src,
                                       const FSClippingFace& planeOwner,
                                       FS_floatT tol = 1e-12);

/* -----------------------------------------------------------------
   ProjectPoly2DTo3D
   ----------------------------------------------------------------- */
std::vector<FSVec3>
ProjectPoly2DTo3D(const FSClippingFace& face,
                  const std::vector<FSVec2>& poly2D);

} // namespace FSClippingUtil

_FS_END_NAMESPACE
#endif // FSCLIPPINGUTIL_H