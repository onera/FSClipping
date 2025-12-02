#ifndef FSCLIPPINGFACE_HPP
#define FSCLIPPINGFACE_HPP

#include "FSBoundingBoxFace.h"
#include "FSClac.h"
#include "FSCommon.h"
#include "FSConfig.h"
#include "FSFace.h"
#include "FSMeshFaceExtractor.h"
#include <FSArray.h>
#include <FSTypes.h>

#include <FSVec3.h>
#include <cassert>
#include <vector>

_FS_BEGIN_NAMESPACE

/**
 * @class FSClippingFace
 * @brief Holds the geometric description of a face that will be used for
 * clipping.
 *
 */
class FSClippingFace {
public:
  /* --------------------------------------------------------------
     Constructors
     -------------------------------------------------------------- */
  explicit FSClippingFace(
      const FSFace &face,
      const FSMeshFaceExtractor
          &faceExtractor) // should we keep face Extractor in the class ?
                          // should we use outFaceNodeCoordinates instead of
                          // faceExtractor ?
      : topo_(face) {
    buildGeometry(faceExtractor);
  }

  /** Build the object from a raw coordinate array (3‑D only). */
  explicit FSClippingFace(const FSFloatArrayT &coords) : topo_() {
    setVerticesFromArray(coords); // triger computeGeometry()
  }

  /** Build the object from a raw std::vector of FSVec3D (3‑D only). */
  explicit FSClippingFace(const std::vector<FSVec3> &vec) : topo_() {
    setVertices(vec); // triger computeGeometry()
  }

  /* --------------------------------------------------------------
     Public read‑only accessors (noexcept because they never throw)
     -------------------------------------------------------------- */
  const FSFace &topo() const noexcept { return topo_; }
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
  /**
   * @brief Build all geometric data (normals, tangents, bounding box,
   *        2‑D projection) from a mesh face.
   *
   * @param face          Original mesh face (only its index/ID is used).
   * @param faceExtractor Helper object that can retrieve node coordinates.
   */
  void buildGeometry(const FSMeshFaceExtractor &faceExtractor);

private:
  /* --------------------------------------------------------------
     Private data members
     -------------------------------------------------------------- */
  const FSFace topo_; // immutable wrapper to the FSDM connectivity
  FSBoundingBoxFace boundingBox_;
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

#endif // FSCLIPPINGFACE_HPP