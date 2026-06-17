#ifndef FSBOUNDARYFACEPROVIDER_H
#define FSBOUNDARYFACEPROVIDER_H

#include <FSMesh.h>

#include "FSClipping/FSClippingFace.h"

#include <unordered_set>

_FS_BEGIN_NAMESPACE

struct QuantizedPoint {
  int x, y, z;

  bool operator==(const QuantizedPoint& other) const { return x == other.x && y == other.y && z == other.z; }

  bool operator<(const QuantizedPoint& other) const
  {
    if(x != other.x)
      return x < other.x;
    if(y != other.y)
      return y < other.y;
    return z < other.z;
  }
};

inline QuantizedPoint Quantize(const FSVec3& p, double tol)
{
  return {static_cast<int>(std::round(p[0] / tol)), static_cast<int>(std::round(p[1] / tol)),
          static_cast<int>(std::round(p[2] / tol))};
}

struct GeomFaceKey {
  std::vector<QuantizedPoint> pts;

  GeomFaceKey(const FSFloatArrayT& coords, double tol)
  {
    const int n = coords.Size(0); // vertex count

    pts.reserve(n);
    for(int i = 0; i < n; ++i) {
      FSVec3 p(coords(i, 0), coords(i, 1), coords(i, 2));
      pts.push_back(Quantize(p, tol));
    }

    Canonicalize();
  }

  void Canonicalize()
  {
    const size_t n = pts.size();
    if(n == 0)
      return;

    std::vector<std::vector<QuantizedPoint> > candidates;

    // forward
    for(size_t i = 0; i < n; ++i) {
      std::vector<QuantizedPoint> rot(n);
      for(size_t j = 0; j < n; ++j)
        rot[j] = pts[(i + j) % n];
      candidates.push_back(rot);
    }

    // reverse
    std::vector<QuantizedPoint> rev(pts.rbegin(), pts.rend());
    for(size_t i = 0; i < n; ++i) {
      std::vector<QuantizedPoint> rot(n);
      for(size_t j = 0; j < n; ++j)
        rot[j] = rev[(i + j) % n];
      candidates.push_back(rot);
    }

    pts = *std::min_element(candidates.begin(), candidates.end());
  }

  bool operator==(const GeomFaceKey& other) const { return pts == other.pts; }
};

struct GeomFaceKeyHash {
  size_t operator()(const GeomFaceKey& k) const
  {
    size_t h = 0;
    for(const auto& p : k.pts) {
      size_t ph = std::hash<int>()(p.x) ^ (std::hash<int>()(p.y) << 1) ^ (std::hash<int>()(p.z) << 2);

      h ^= ph + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};

struct BoundaryExtraction {
  std::vector<FSClippingFace> faces;
  std::unordered_set<GeomFaceKey, GeomFaceKeyHash> faceKeys;

  std::unordered_map<FS_intT, std::set<FS_intT> > volumeCells;
  std::unordered_map<FS_intT, std::set<FS_intT> > surfaceCells;
};

class FSBoundaryFaceProvider
{
public:
  static BoundaryExtraction Extract(FSMesh& mesh, FSMeshFaceExtractor& ex, FS_intT boundaryMarker, FS_floatT tol,
                                    bool matchRemoteFaces = true);
};


_FS_END_NAMESPACE

#endif