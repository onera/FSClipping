#ifndef FSCLIPPINGFACE_H
#define FSCLIPPINGFACE_H

#include "FSFace.h"
#include "FSMeshFaceExtractor.h"

#include <atomic>

_FS_BEGIN_NAMESPACE

/* --------------------------------------------------------------
   Bounding box of a face
   -------------------------------------------------------------- */
struct FSBoundingBoxFace {
  FS_floatT boxMinMax[6];

  FSBoundingBoxFace()
    : boxMinMax(+FS_FLOATT_MAX, +FS_FLOATT_MAX, +FS_FLOATT_MAX, -FS_FLOATT_MAX, -FS_FLOATT_MAX, -FS_FLOATT_MAX)
  {
  }

  void ExpandToInclude(const FSVec3& p)
  {
    boxMinMax[0] = FSMin(boxMinMax[0], p[0]);
    boxMinMax[1] = FSMin(boxMinMax[1], p[1]);
    boxMinMax[2] = FSMin(boxMinMax[2], p[2]);
    boxMinMax[3] = FSMax(boxMinMax[3], p[0]);
    boxMinMax[4] = FSMax(boxMinMax[4], p[1]);
    boxMinMax[5] = FSMax(boxMinMax[5], p[2]);
  }
};

class FSClippingFace
{
public:
  /* --------------------------------------------------------------
     Constructors
     -------------------------------------------------------------- */
  explicit FSClippingFace(const FSFace& face, const FS_intT faceIndex, const FSFloatArrayT& faceNodeCoordinates)
    : topo_(face), faceIndex_(faceIndex)
  {
    buildGeometry(faceNodeCoordinates);
  };

  explicit FSClippingFace(const FSFloatArrayT& coords) : topo_(), faceIndex_(++s_nextId_)
  {
    setVerticesFromArray(coords); // triggers computeGeometry()
  }

  explicit FSClippingFace(const std::vector<FSVec3>& vec) : topo_(), faceIndex_(++s_nextId_)
  {
    setVertices(vec); // triggers computeGeometry()
  }

  /* --------------------------------------------------------------
     Public read-only accessors
     -------------------------------------------------------------- */
  const FSFace& topo() const noexcept { return topo_; }
  const FS_intT& faceIndex() const noexcept { return faceIndex_; }

  // Reset the auto-increment id used by the coordinate-only constructors. Only
  // meant for unit tests: s_nextId_ is a process-wide static that keeps counting
  // across TEST bodies run in the same process, so a test asserting on absolute
  // faceIndex values must reset it in its fixture SetUp to be independent of run
  // order / execution mode.
  static void ResetIdCounter() noexcept { s_nextId_ = 0; }
  const FSBoundingBoxFace& boundingBox() const noexcept { return boundingBox_; }
  const std::vector<FSVec3>& vertices() const noexcept { return vertices_; }
  const std::vector<FSVec2>& projected2D() const noexcept { return projected2D_; }
  const FSVec3& normal() const noexcept { return normal_; }
  const FSVec3& tangentVector1() const noexcept { return tangentVector1_; }
  const FSVec3& tangentVector2() const noexcept { return tangentVector2_; }

  /* --------------------------------------------------------------
     Public mutators
     -------------------------------------------------------------- */
  void setVertices(const std::vector<FSVec3>& verts)
  {
    vertices_ = verts;
    computeGeometry();
  }

  void setVerticesFromArray(const FSFloatArrayT& coords)
  {
    const FS_intT size = coords.Size() / FS_3D;
    vertices_.clear();
    vertices_.reserve(static_cast<std::size_t>(size));

    for(FS_intT i = 0; i < size; ++i)
      vertices_.emplace_back(coords(i, 0), coords(i, 1), coords(i, 2));
    computeGeometry();
  }

  /* --------------------------------------------------------------
     Full geometry construction
     -------------------------------------------------------------- */
  void buildGeometry(const FSFloatArrayT& faceNodeCoordinates);

private:
  /* --------------------------------------------------------------
     Data members
     -------------------------------------------------------------- */
  const FSFace topo_;                    // immutable wrapper to the FSDM connectivity
  const FS_intT faceIndex_;              // global index of the face
  FSBoundingBoxFace boundingBox_;        // axis-aligned bounding box
  static std::atomic<FS_intT> s_nextId_; // auto-increment id used by coordinate-only constructors (unit tests)
  std::vector<FSVec3> vertices_;         // 3D vertices
  std::vector<FSVec2> projected2D_;      // 2D projection onto the local basis
  FSVec3 normal_;                        // unit normal
  FSVec3 tangentVector1_;                // first orthonormal tangent
  FSVec3 tangentVector2_;                // second orthonormal tangent

  /* Compute normal, orthonormal basis, bounding box and 2D projection.
     Called whenever the vertex list changes. */
  void computeGeometry();
};

_FS_END_NAMESPACE

#endif // FSCLIPPINGFACE_H
