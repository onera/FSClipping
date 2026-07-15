#ifndef FSFACEMATCHER_H
#define FSFACEMATCHER_H

#include "FSBVHTree.h"
#include "FSBoundingBoxUtil.h"

#include "FSClipping/FSClippingFace.h"

_FS_BEGIN_NAMESPACE

struct FSFaceMatch {

  // --- Identification of the two faces (local indices in their arrays) ---
  FS_intT face1 = -1; // index in the FSClippingFace list of mesh 1
  FS_intT face2 = -1; // index in the FSClippingFace list of mesh 2 (-1 for NOT_COVERED matches)

  // --- Ownership information (already known from FSDM / MPI) ---
  FS_intT elemOwner1 = -1; // element owner of the face 1
  FS_intT elemOwner2 = -1; // element owner of the face 2

  FS_intT faceOwner1 = -1; // face owner of the face 1
  FS_intT faceOwner2 = -1; // face owner of the face 2

  FSMeshEnums::CellType elemOwnerType1 = FSMeshEnums::CT_Undefined;
  FSMeshEnums::CellType elemOwnerType2 = FSMeshEnums::CT_Undefined;

  FSMeshEnums::CellType faceOwnerType1 = FSMeshEnums::CT_Undefined;
  FSMeshEnums::CellType faceOwnerType2 = FSMeshEnums::CT_Undefined;

  // --- Type of geometric relation ---
  // NOT_COVERED: the subject face is not entirely covered by the clipped mesh
  // (e.g. rim faces of two cylinders in relative rotation). The face is kept
  // unclipped: the match carries the original face polygon and has no face2.
  enum MatchType : FS_intT {
    UNKNOWN = 0,
    IDENTICAL = 1,
    INCLUDED = 2,
    INTERSECTING = 3,
    NOT_COVERED = 4
  } type = UNKNOWN;

  // --- Geometric measure ---
  FS_floatT intersectedArea = 0.0;

  // --- clippedPoly2D polygon in the local plane of face1 ---
  // stored in 2D because polygon clipping is performed in 2D
  std::vector<FSVec2> clippedPoly2D;

  // Intersection polygon re-projected in 3D in the frame of face1
  std::vector<FSVec3> clippedPoly3D;

  // Intersection polygon re-projected in 3D in the frame of face2.
  // ComputeInvertedMatches swaps the two so that clippedPoly3D always refers
  // to the current face1 frame. This ensures consistent float values when the
  // same geometric vertex appears in multiple matches for the same cell.
  std::vector<FSVec3> clippedPoly3D_face2;

  // --- Global numbering (level 2 parallelism) ---
  // Assigned by the clipper proc (AssignGlobalNodeIds) before scattering the
  // matches, so that two mesh procs sharing a node at their partition
  // interface receive the same global node ID by construction.
  std::vector<FS_intT> nodeGlobalIds; // parallel to clippedPoly3D
  FS_intT globalCellId = -1;          // global ID of the future Poly2D cell

  FSClac::sizeT GetBufSize(FSClac& clac) const;
  void Pack(FSClac& clac);
  void Unpack(FSClac& clac);
};

class FSFaceMatcher
{
public:
  FSFaceMatcher(FSClac& clacClipped, const std::vector<FSClippingFace>& subjectFaces,
                const std::vector<FSClippingFace>& clippedFaces, FS_floatT tol = 1e-12)
    : subjectFaces_(subjectFaces), clippedFaces_(clippedFaces), tol_(tol)
  {

    FS_sizeT size = static_cast<FS_sizeT>(clippedFaces.size());
    FSFloatArrayT boundingBoxes(size, 6);
    FSIntArrayT indices(size);
    for(FS_sizeT i = 0; i < size; ++i) {
      for(FS_intT j = 0; j < 6; j++)
        boundingBoxes(i, j) = clippedFaces[i].boundingBox().boxMinMax[j];
      indices(i) = clippedFaces[i].faceIndex(); // global ID
    }
    if(!FSBoundingBoxUtil::GatherBoundingBoxesIntoBVHTree(&clacClipped, boundingBoxes, indices, bvhClipped_))
      FSError.SetAndPrintAndExit("Failed to build BVH for clipped face");
  }

  // Computes the match list for the subject side: raw face-pair intersections,
  // then subject faces not fully covered by the clipped mesh are kept whole as
  // single NOT_COVERED matches (see PreserveUncoveredFaces).
  void ComputeMatches(std::vector<FSFaceMatch>& outMatches);

  // Computes the match list for the clipped side (face1/face2 roles swapped).
  // Rebuilt from the raw matches — not from the subject-side list — so that the
  // coverage decision is made independently on each side: a partial clip
  // dropped for an uncovered subject face is still a valid intersection for
  // the clipped face it belongs to. Requires ComputeMatches to have run.
  void ComputeInvertedMatches(std::vector<FSFaceMatch>& outMatches) const;

  // Assigns a globally consistent numbering across one mesh side:
  //  - nodeGlobalIds: geometric deduplication (NodeKey, tol) of every point of
  //    every clippedPoly3D — two points within tol get the same ID, IDs are
  //    contiguous 0..K-1;
  //  - globalCellId: match index 0..n-1 (global ID of the future Poly2D cell).
  // Called by the clipper on each side's final match list before scattering,
  // so procs sharing an interface node receive the same ID by construction.
  // Returns the number of unique node IDs assigned (K).
  static FS_intT AssignGlobalNodeIds(std::vector<FSFaceMatch>& matches, FS_floatT tol);

private:
  const std::vector<FSClippingFace>& subjectFaces_;
  const std::vector<FSClippingFace>& clippedFaces_;

  FSBVHTree bvhClipped_;
  FS_floatT tol_;

  // All face-pair intersections, before any coverage filtering. Kept so that
  // ComputeInvertedMatches can apply the coverage criterion on the clipped side.
  std::vector<FSFaceMatch> rawMatches_;

  bool ComputeMatch(const FSClippingFace& f1, const FSClippingFace& f2, FSFaceMatch& out) const;

  // For every face of `faces` whose matches (identified via face1) do not sum
  // up to the full face area, drops the partial clips and appends one
  // NOT_COVERED match carrying the original face polygon instead.
  void PreserveUncoveredFaces(std::vector<FSFaceMatch>& matches, const std::vector<FSClippingFace>& faces) const;
};

_FS_END_NAMESPACE
#endif // FSFACEMATCHER_H
