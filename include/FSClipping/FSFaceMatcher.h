#ifndef FSFACEMATCHER_H
#define FSFACEMATCHER_H

#include "FSBVHTree.h"
#include "FSBoundingBoxUtil.h"

#include "FSClipping/FSClippingFace.h"

_FS_BEGIN_NAMESPACE

struct FSFaceMatch {

  // --- Identification of the two faces (local indices in their arrays) ---
  FS_intT face1; // index in the FSClippingFace list of mesh 1
  FS_intT face2; // index in the FSClippingFace list of mesh 2

  // --- Ownership information (already known from FSDM / MPI) ---
  FS_intT elemOwner1; // element owner of the face 1
  FS_intT elemOwner2; // element owner of the face 2

  FSMeshEnums::CellType elemOwnerType1;
  FSMeshEnums::CellType elemOwnerType2;

  // --- Type of geometric relation ---
  enum MatchType : FS_intT {
    UNKNOWN = 0,
    IDENTICAL = 1,
    INCLUDED = 2,
    INTERSECTING = 3
  } type = UNKNOWN;

  // --- Geometric measure ---
  // double FS_IntTersectedArea = 0.0;

  // --- clippedPoly2D polygon in the local plane of face1 ---
  // stored in 2D because polygon clipping is performed in 2D
  std::vector<FSVec2> clippedPoly2D;

  // FS_IntTersection polygon re‑projected in 3D in face 1
  std::vector<FSVec3> clippedPoly3D;
};

class FSFaceMatcher
{ // does we need the faceExtractor for the parralel
  // computation ?
public:
  FSFaceMatcher(
    FSClac& clacClipped, // const FSMeshFaceExtractor &faceExtractor,
    const std::vector<FSClippingFace>& subjectFaces,
    const std::vector<FSClippingFace>& clippedFaces)
    : subjectFaces_(subjectFaces),
      clippedFaces_(clippedFaces)
  { //, faceExtractorClipped_(faceExtractor) {

    FS_sizeT size = static_cast<FS_sizeT>(clippedFaces.size());
    FSFloatArrayT boundingBoxes(size, 6);
    FSIntArrayT indices(size);
    for(FS_sizeT i = 0; i < size; ++i) {
      for(FS_intT j = 0; j < 6; j++)
        boundingBoxes(i, j) = clippedFaces[i].boundingBox().boxMinMax[j];
      indices(i) = clippedFaces[i].faceIndex(); // global ID
    }
    if(!FSBoundingBoxUtil::GatherBoundingBoxesIntoBVHTree(
         &clacClipped, boundingBoxes, indices, bvhClipped_))
      FSError.SetAndPrintAndExit("Failed to build BVH for clipped face");
  }

  void ComputeMatches(std::vector<FSFaceMatch>& outMatches);

private:
  const std::vector<FSClippingFace>& subjectFaces_;
  const std::vector<FSClippingFace>& clippedFaces_;
  //  const FSMeshFaceExtractor &faceExtractorClipped_;
  FSBVHTree bvhClipped_;

  bool ComputeMatch(const FSClippingFace& f1, const FSClippingFace& f2,
                    FSFaceMatch& out) const;
};

_FS_END_NAMESPACE
#endif // FSFACEMATCHER_H
