#ifndef FSFACEMATCHER_H
#define FSFACEMATCHER_H

#include "FSBVHTree.h"
#include "FSBoundingBoxUtil.h"

#include "FSClipping/FSClippingFace.h"

_FS_BEGIN_NAMESPACE

struct FSFaceMatch {

  // Local index of the two matches faces
  FS_intT face1 = -1;
  FS_intT face2 = -1;

  // Local index of the element owner of the two faces
  FS_intT elemOwner1 = -1;
  FS_intT elemOwner2 = -1;

  // ProcID of each faces (local to their own clac)
  FS_intT ownerProc1 = -1;
  FS_intT ownerProc2 = -1;

  FS_intT faceOwner1 = -1;
  FS_intT faceOwner2 = -1;

  FSMeshEnums::CellType elemOwnerType1 = FSMeshEnums::CT_Undefined;
  FSMeshEnums::CellType elemOwnerType2 = FSMeshEnums::CT_Undefined;

  FSMeshEnums::CellType faceOwnerType1 = FSMeshEnums::CT_Undefined;
  FSMeshEnums::CellType faceOwnerType2 = FSMeshEnums::CT_Undefined;

  // --- Type of geometric relation ---
  enum MatchType : FS_intT {
    UNKNOWN = 0,
    IDENTICAL = 1,
    INCLUDED = 2,
    INTERSECTING = 3,
    NOT_COVERED = 4 // could be if the match of face (A/B) are note totaly recover by the mesh (A/B)
  } type = UNKNOWN;

  FS_floatT intersectedArea = 0.0;

  // 2D polygon resulted by the match in the local plane of face1
  std::vector<FSVec2> clippedPoly2D;

  // 3D polygon re-projected in the frame of face1
  std::vector<FSVec3> clippedPoly3D;

  // 3D polygon re-projected in the frame of face2.
  std::vector<FSVec3> clippedPoly3D_face2;

  // Global numbering of each node of a match
  std::vector<FS_intT> nodeGlobalIds;
  // On peut pas creer cette liste dans le Cell2NodeBuilder ?

  // Global ID of the future Poly2D cell
  FS_intT globalCellId = -1;
  // Meme remarque, la numerotation globale ne peut pas etre faite dans Cell2NodeBuilder ?

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


  void ComputeMatches(std::vector<FSFaceMatch>& outMatches);

  void ComputeInvertedMatches(std::vector<FSFaceMatch>& outMatches) const;

  /*
  Assigns a globally consistent numbering across one mesh side:
    - nodeGlobalIds: geometric deduplication (NodeKey, tol) of every point of every clippedPoly3D
    - globalCellId: match index 0..n-1 (global ID of the future Poly2D cell)
  */
  static FS_intT AssignGlobalNodeIds(std::vector<FSFaceMatch>& matches, FS_floatT tol);

private:
  const std::vector<FSClippingFace>& subjectFaces_;
  const std::vector<FSClippingFace>& clippedFaces_;

  FSBVHTree bvhClipped_;
  FS_floatT tol_;

  // All face-pair intersections, before any coverage filtering.
  std::vector<FSFaceMatch> rawMatches_;

  bool ComputeMatch(const FSClippingFace& f1, const FSClippingFace& f2, FSFaceMatch& out) const;

  // Preserve the faces of the clipping algorithm when they are not totaly covered by the facing mesh
  void PreserveUncoveredFaces(std::vector<FSFaceMatch>& matches, const std::vector<FSClippingFace>& faces) const;
};

_FS_END_NAMESPACE
#endif // FSFACEMATCHER_H
