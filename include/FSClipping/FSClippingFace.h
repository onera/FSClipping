#ifndef FSCLIPPINGFACE_H
#define FSCLIPPINGFACE_H

#include "FSFace.h"
#include "FSMeshFaceExtractor.h"

#include <atomic>

_FS_BEGIN_NAMESPACE

struct FSBoundingBoxFace {
  FS_floatT boxMinMax[6];

  FSBoundingBoxFace()
      : boxMinMax(+FS_FLOATT_MAX, +FS_FLOATT_MAX, +FS_FLOATT_MAX,
                  -FS_FLOATT_MAX, -FS_FLOATT_MAX, -FS_FLOATT_MAX) {}

  void ExpandToInclude(const FSVec3 &p) {
    boxMinMax[0] = FSMin(boxMinMax[0], p[0]);
    boxMinMax[1] = FSMin(boxMinMax[1], p[1]);
    boxMinMax[2] = FSMin(boxMinMax[2], p[2]);
    boxMinMax[3] = FSMax(boxMinMax[3], p[0]);
    boxMinMax[4] = FSMax(boxMinMax[4], p[1]);
    boxMinMax[5] = FSMax(boxMinMax[5], p[2]);
  }
};

class FSClippingFace {
public:
  /* --------------------------------------------------------------
     Constructors
     -------------------------------------------------------------- */
  explicit FSClippingFace(const FSFace &face, const FS_intT faceIndex,
                          const FSFloatArrayT &faceNodeCoordinates)
      : topo_(face), faceIndex_(faceIndex) {
    buildGeometry(faceNodeCoordinates);
  };

  /** Build the object from a raw coordinate array (3‑D only). */
  explicit FSClippingFace(const FSFloatArrayT &coords)
      : topo_(), faceIndex_(++s_nextId_) {
    setVerticesFromArray(coords); // triger computeGeometry()
  }

  /** Build the object from a raw std::vector of FSVec3D (3‑D only). */
  explicit FSClippingFace(const std::vector<FSVec3> &vec)
      : topo_(), faceIndex_(++s_nextId_) {
    setVertices(vec); // triger computeGeometry()
  }

  /* --------------------------------------------------------------
     Public read‑only accessors (noexcept because they never throw)
     -------------------------------------------------------------- */
  const FSFace &topo() const noexcept { return topo_; }
  const FS_intT &faceIndex() const noexcept { return faceIndex_; }
  const FSBoundingBoxFace &boundingBox() const noexcept { return boundingBox_; }
  const std::vector<FSVec3> &vertices() const noexcept { return vertices_; }
  const std::vector<FSVec2> &projected2D() const noexcept {
    return projected2D_;
  }
  const FSVec3 &normal() const noexcept { return normal_; }
  const FSVec3 &tangentVector1() const noexcept { return tangentVector1_; }
  const FSVec3 &tangentVector2() const noexcept { return tangentVector2_; }

  /* --------------------------------------------------------------
     Public mutators – thin wrappers that keep the object in a valid state
     -------------------------------------------------------------- */

  /** Replace the whole vertex list (geometry will be recomputed). */
  void setVertices(const std::vector<FSVec3> &verts) {
    vertices_ = verts;
    // recompute everything from the new vertices
    computeGeometry();
  }

  /** Fill the vertex list from a raw FSFloatArrayT (expects 3‑D data). */
  void setVerticesFromArray(const FSFloatArrayT &coords) {
    const FS_intT size = coords.Size() / FS_3D;
    vertices_.clear();
    vertices_.reserve(static_cast<std::size_t>(size));

    for (FS_intT i = 0; i < size; ++i)
      vertices_.emplace_back(coords(i, 0), coords(i, 1), coords(i, 2));
    // recompute geometry now that we have the vertices
    computeGeometry();
  }

  /* --------------------------------------------------------------
     Public wrapper that performs the full geometry construction
     -------------------------------------------------------------- */
  void buildGeometry(const FSFloatArrayT &faceNodeCoordinates);

private:
  /* --------------------------------------------------------------
     Private data members
     -------------------------------------------------------------- */
  const FSFace topo_; // immutable wrapper to the FSDM connectivity, should I
                      // keep the topology ?
  const FS_intT faceIndex_;
  FSBoundingBoxFace boundingBox_;
  static std::atomic<FS_intT> s_nextId_;
  std::vector<FSVec3> vertices_;    // 3‑D vertices of the face
  std::vector<FSVec2> projected2D_; // 2‑D projection in the local basis
  FSVec3 normal_;                   // unit normal
  FSVec3 tangentVector1_;           // first orthonormal basis vector
  FSVec3 tangentVector2_;           // second orthonormal basis vector

  /* --------------------------------------------------------------
     Private implementation details
     -------------------------------------------------------------- */
  /** Compute normal, orthonormal basis, bounding box and 2‑D projection.
      Called internally whenever the vertex list changes. */
  void computeGeometry();
};

_FS_END_NAMESPACE

#endif // FSCLIPPINGFACE_H