#include "FSClipping/FSClippingUtil.h"
#include "FSClipping/FSClippingFace.h"

_FS_BEGIN_NAMESPACE
namespace FSClippingUtil {

bool AreFacesIdentical(const FSClippingFace& f1,
                       const FSClippingFace& f2,
                       FS_floatT tol)
{
  // improve this algorithm to have nlog(n) cost

  const auto& v1 = f1.vertices();
  const auto& v2 = f2.vertices();

  if(v1.size() != v2.size())
    return false;

  const auto n = v1.size();
  std::vector<bool> used(n, false);

  for(std::size_t i = 0; i < n; ++i) {
    const auto& p = v1[i];
    bool matched = false;
    for(std::size_t j = 0; j < n; ++j) {
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

bool ClipConvexPolygon(const std::vector<FSVec2>& subject,
                       const std::vector<FSVec2>& clipper,
                       std::vector<FSVec2>& output,
                       FS_floatT tol)
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

std::vector<FSVec2> ProjectFaceInPlane(const FSClippingFace& src,
                                       const FSClippingFace& planeOwner,
                                       FS_floatT tol)
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

std::vector<FSVec3>
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