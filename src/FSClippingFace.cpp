#include "FSClipping/FSClippingFace.h"

_FS_BEGIN_NAMESPACE

std::atomic<FS_intT> FSClippingFace::s_nextId_{0};
/* -----------------------------------------------------------------
   Implementation of the public wrapper `buildGeometry`.
   This function simply forwards the call to the private
   `computeGeometry` after having retrieved the vertex coordinates
   from the mesh extractor.
   ----------------------------------------------------------------- */
void FSClippingFace::buildGeometry(const FSFloatArrayT& faceNodeCoordinates)
{
  // -----------------------------------------------------------------
  // 1) Retrieve coordinates of the face nodes (identical to the legacy code)
  // -----------------------------------------------------------------
  assert(topo_._faceFSDM);

  //  FSFloatArrayT outFaceNodeCoordinates;
  //  const bool orientation =
  //      faceExtractor.GetFaceNodeCoordinates(faceIndex_,
  //      outFaceNodeCoordinates);
  //
  //  const FS_intT size = outFaceNodeCoordinates.Size() / FS_3D;
  //  vertices_.clear();
  //  vertices_.reserve(static_cast<std::size_t>(size));
  //
  const FS_intT size = faceNodeCoordinates.Size() / FS_3D;
  for(FS_intT i = 0; i < size; ++i)
    vertices_.emplace_back(faceNodeCoordinates(i, 0), faceNodeCoordinates(i, 1),
                           faceNodeCoordinates(i, 2));

  // -----------------------------------------------------------------
  // 2) Compute geometry (normal, tangents, bounding box, projection)
  // -----------------------------------------------------------------
  computeGeometry();

  // Apply orientation correction if the extractor reported a reversed order
  // if (!orientation) {
  //   normal_ = -normal_;
  //   // Re‑compute the second tangent because it depends on the normal
  //   tangentVector2_ = normal_.CrossProduct(tangentVector1_);
  // }
}

/* -----------------------------------------------------------------
   Private helper that does the heavy‑lifting.  It assumes that
   `vertices_` already contains the 3‑D points of the face.
   ----------------------------------------------------------------- */
void FSClippingFace::computeGeometry()
{
  // -------------------------------------------------------------
  // 1) Tangent 1  (v1 - v0)  → normalized
  // -------------------------------------------------------------
  FSVec3 t1 = vertices_[1] - vertices_[0];
  tangentVector1_ = t1 / t1.L2Norm();

  // -------------------------------------------------------------
  // 2) Normal = (t1 × (v2 - v0)).normalized()
  // -------------------------------------------------------------
  FSVec3 t2 = vertices_[2] - vertices_[0];
  normal_ = t1.CrossProduct(t2);
  normal_ /= normal_.L2Norm();

  // -------------------------------------------------------------
  // 3) Tangent 2 = normal × tangentVector1 (already orthonormal)
  // -------------------------------------------------------------
  tangentVector2_ = normal_.CrossProduct(tangentVector1_);

  // -------------------------------------------------------------
  // 4) Bounding box
  // -------------------------------------------------------------
  boundingBox_ = {}; // reset
  for(const auto& v : vertices_)
    boundingBox_.ExpandToInclude(v);

  // -------------------------------------------------------------
  // 5) 2‑D projection of each vertex on the local basis
  // -------------------------------------------------------------
  projected2D_.clear();
  projected2D_.reserve(vertices_.size());

  const FSVec3 origin = vertices_[0];
  for(const FSVec3& P : vertices_) {
    FSVec3 d = P - origin;
    FS_floatT u = d.InnerProduct(tangentVector1_);
    FS_floatT v = d.InnerProduct(tangentVector2_);
    projected2D_.emplace_back(u, v);
  }
}

_FS_END_NAMESPACE
