#ifndef FSCLIPPINGUTIL_H
#define FSCLIPPINGUTIL_H

#include "FSPlaneSurface.h"

#include "FSClipping/FSClippingFace.h"

_FS_BEGIN_NAMESPACE

namespace FSClippingUtil {

/* -----------------------------------------------------------------
   AreFacesIdentical
   ----------------------------------------------------------------- */
inline bool AreFacesIdentical(const FSClippingFace& f1,
                              const FSClippingFace& f2,
                              FS_floatT tol = static_cast<FS_floatT>(1e-12))
{
  // improve this algorithm to have nlog(n) cost

  const auto& v1 = f1.vertices();
  const auto& v2 = f2.vertices();

  if(v1.size() != v2.size())
    return false;

  const FS_intT n = static_cast<FS_intT>(v1.size());
  std::vector<bool> used(n, false);

  for(FS_intT i = 0; i < n; ++i) {
    const auto& p = v1[i];
    bool matched = false;
    for(FS_intT j = 0; j < n; ++j) {
      if(used[j])
        continue;
      if((p - v2[j]).L2Norm() < tol) {
        used[j] = true;
        matched = true;
        break;
      }
    }
    if(!matched)
      return false;
  }
  return true;
}

/* -----------------------------------------------------------------
   IsPointInsideConvexPolygon
   ----------------------------------------------------------------- */
inline bool
IsPointInsideConvexPolygon(const FSVec2& P, const std::vector<FSVec2>& poly,
                           FS_floatT tol = static_cast<FS_floatT>(1e-12))
{
  const FS_intT n = static_cast<FS_intT>(poly.size());
  for(FS_intT i = 0; i < n; ++i) {
    const FSVec2& A = poly[i];
    const FSVec2& B = poly[(i + 1) % n];
    FSVec2 edge = B - A;
    // normale intérieure (orientation CCW)
    FSVec2 normal(-edge[1], edge[0]);
    FS_floatT d = (P - A).InnerProduct(normal);
    if(d < -tol)
      return false; // à l’extérieur
  }
  return true;
}

/* -----------------------------------------------------------------
   IsFaceInsideTheOther
   ----------------------------------------------------------------- */
inline bool
IsFaceInsideTheOther(const std::vector<FSVec2>& f1,
                     const std::vector<FSVec2>& f2,
                     FS_floatT tol = static_cast<FS_floatT>(1e-12))
{
  for(const FSVec2& P : f1) {
    if(!IsPointInsideConvexPolygon(P, f2, tol))
      return false;
  }
  return true;
}

/* -----------------------------------------------------------------
   IsPointInsideEdge
   ----------------------------------------------------------------- */
inline bool IsPointInsideEdge(const FSVec2& P, const FSVec2& A, const FSVec2& B,
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
inline FSVec2 ComputeIntersection(const FSVec2& P1, const FSVec2& P2,
                                  const FSVec2& A, const FSVec2& B)
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
   ClipConvexPolygon
   ----------------------------------------------------------------- */
inline bool ClipConvexPolygon(const std::vector<FSVec2>& subject,
                              const std::vector<FSVec2>& clipper,
                              std::vector<FSVec2>& output,
                              FS_floatT tol = static_cast<FS_floatT>(1e-12))
{
  // Validate input polygons
  if(subject.size() < 3 || clipper.size() < 3)
    return false;

  output = subject;
  std::vector<FSVec2> input;

  // Iterate over each edge of the clipping polygon
  for(size_t i = 0; i < clipper.size(); ++i) {
    const FSVec2& A = clipper[i];
    const FSVec2& B = clipper[(i + 1) % clipper.size()];

    input = std::move(output);
    output.clear();
    output.reserve(input.size()); // reserve to avoid reallocations

    if(input.empty())
      return false;

    // We control the push back by checking if the new vertice is not the same
    // than the last
    auto push_unique = [&](const FSVec2& pt) {
      if(output.empty() || (pt - output.back()).Norm() > tol * tol)
        output.push_back(pt);
    };

    // Process each edge of the current subject polygon
    for(size_t j = 0; j < input.size(); ++j) {
      const FSVec2& P1 = input[j];
      const FSVec2& P2 = input[(j + 1) % input.size()];

      bool P1_inside = IsPointInsideEdge(P1, A, B, tol);
      bool P2_inside = IsPointInsideEdge(P2, A, B, tol);

      if(P1_inside && P2_inside) {
        // Both vertices inside – keep the end vertex
        push_unique(P2);
      } else if(P1_inside && !P2_inside) {
        // Edge leaves the clipping region – keep intersection only
        push_unique(ComputeIntersection(P1, P2, A, B));
      } else if(!P1_inside && P2_inside) {
        // Edge enters the clipping region – keep intersection and end vertex
        push_unique(ComputeIntersection(P1, P2, A, B));
        push_unique(P2);
      }
      // else both outside – nothing to add
    }
  }
  return output.size() >= 3;
}

/* -----------------------------------------------------------------
   ProjectFaceInPlane
   ----------------------------------------------------------------- */
inline std::vector<FSVec2> ProjectFaceInPlane(const FSClippingFace& src,
                                              const FSClippingFace& planeOwner,
                                              FS_floatT tol = 1e-12)
{
  FSPlaneSurface pl;
  pl.Init(planeOwner.vertices()[0], planeOwner.tangentVector1(),
          planeOwner.tangentVector2());

  std::vector<FSVec2> proj;
  proj.reserve(src.vertices().size());

  for(const FSVec3& P : src.vertices()) {
    FSVec3 dummy3;
    FSVec2 uv;
    pl.Project(P, dummy3, uv, tol);
    proj.emplace_back(uv);
  }
  return proj;
}

/* -----------------------------------------------------------------
   ProjectPoly2DTo3D
   ----------------------------------------------------------------- */
inline std::vector<FSVec3>
ProjectPoly2DTo3D(const FSClippingFace& face,
                  const std::vector<FSVec2>& poly2D)
{
  std::vector<FSVec3> poly3D;
  poly3D.reserve(poly2D.size());

  const FSVec3& O = face.vertices()[0];     // origine du plan
  const FSVec3& e1 = face.tangentVector1(); // axe u
  const FSVec3& e2 = face.tangentVector2(); // axe v

  for(const FSVec2& uv : poly2D) {
    // P = O + u·e1 + v·e2
    poly3D.emplace_back(O + e1 * uv[0] + e2 * uv[1]);
  }
  return poly3D;
}

} // namespace FSClippingUtil

_FS_END_NAMESPACE
#endif // FSCLIPPINGUTIL_H