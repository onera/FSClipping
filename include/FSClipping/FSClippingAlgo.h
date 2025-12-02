#ifndef FSCLIPPINGALGO_HPP
#define FSCLIPPINGALGO_HPP

#include "FSClippingFace.h"
#include "FSClippingUtil.h"

_FS_BEGIN_NAMESPACE

namespace FSClippingAlgo {

static bool IsPointInsideEdge(const FSVec2 &P, const FSVec2 &A, const FSVec2 &B,
                              FS_floatT tol = static_cast<FS_floatT>(1e-12)) {
  // (B-A) × (p-A) >= 0  →  p is on the left of AB
  FS_floatT cross =
      (B[0] - A[0]) * (P[1] - A[1]) - (B[1] - A[1]) * (P[0] - A[0]);
  return cross >= -tol;
}

static FSVec2 ComputeIntersection(const FSVec2 &P1, const FSVec2 &P2,
                                  const FSVec2 &A, const FSVec2 &B) {
  FSVec2 U = P2 - P1; // direction of subject segment
  FSVec2 V = B - A;   // direction of clipping edge

  FS_floatT denom = U[0] * V[1] - U[1] * V[0];
  if (std::abs(denom) < static_cast<FS_floatT>(1e-12))
    return P1; // Parallel – return one endpoint as fallback

  FSVec2 AP = A - P1;
  FS_floatT t = (AP[0] * V[1] - AP[1] * V[0]) / denom;
  return P1 + U * t;
}

static bool ClipConvexPolygon(const std::vector<FSVec2> &subject,
                              const std::vector<FSVec2> &clipper,
                              std::vector<FSVec2> &output,
                              FS_floatT tol = static_cast<FS_floatT>(1e-12)) {
  // Validate input polygons
  if (subject.size() < 3 || clipper.size() < 3)
    return false;

  output = subject;
  std::vector<FSVec2> input;

  // Iterate over each edge of the clipping polygon
  for (size_t i = 0; i < clipper.size(); ++i) {
    const FSVec2 &A = clipper[i];
    const FSVec2 &B = clipper[(i + 1) % clipper.size()];

    input = std::move(output);
    output.clear();
    output.reserve(input.size()); // reserve to avoid reallocations

    if (input.empty())
      return false;

    // Process each edge of the current subject polygon
    for (size_t j = 0; j < input.size(); ++j) {
      const FSVec2 &P1 = input[j];
      const FSVec2 &P2 = input[(j + 1) % input.size()];

      bool P1_inside = IsPointInsideEdge(P1, A, B, tol);
      bool P2_inside = IsPointInsideEdge(P2, A, B, tol);

      if (P1_inside && P2_inside) {
        // Both vertices inside – keep the end vertex
        output.push_back(P2);
      } else if (P1_inside && !P2_inside) {
        // Edge leaves the clipping region – keep intersection only
        output.push_back(ComputeIntersection(P1, P2, A, B));
      } else if (!P1_inside && P2_inside) {
        // Edge enters the clipping region – keep intersection and end vertex
        output.push_back(ComputeIntersection(P1, P2, A, B));
        output.push_back(P2);
      }
      // else both outside – nothing to add
    }
  }
  return output.size() >= 3;
}

} // namespace FSClippingAlgo

_FS_END_NAMESPACE
#endif // FSCLIPPINGUTIL_HPP