#ifndef FSPOLYFACEBUILDER_H
#define FSPOLYFACEBUILDER_H

#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSCell2NodeBuilder.h"

_FS_BEGIN_NAMESPACE

struct FaceData {
  std::vector<FS_intT> nodeIds;
};


class FSPolyFaceBuilder
{
public:
  FSPolyFaceBuilder(const std::vector<FSFaceMatch>& matches,
                    const FSCell2NodeBuilder& cell2NodeBuilder,
                    FS_floatT tol = 1e-12) : matches_(matches),
                                             cell2NodeBuilder_(cell2NodeBuilder),
                                             tol_(tol) {};
  // FSPolyFaceBuilder(std::vector<std::vector<FaceData> >& cellFaces,
  //                   const FSCell2NodeBuilder& cell2NodeBuilder,
  //                   FS_floatT tol = 1e-12) : cellFaces_(cellFaces), cell2NodeBuilder_(cell2NodeBuilder), tol_(tol) {};

  void AddInnerFaces(FSMeshEnums::CellType,
                     const FSCellPool& cellPool,
                     FS_intT cell,
                     FS_intT face,
                     const FSFloatArrayT& oldCoords);

  void Build(FSMeshPolyFaceStorage& polyFaces); // build the FSMeshPolyFaceStorage for FSDM

  void CollectMatchesFaces(); // main function, construct the cellFaces containers with the matches

  const std::vector<std::vector<FaceData> >& CellFaces() const noexcept { return cellFaces_; };

private:
  std::vector<std::vector<FaceData> > cellFaces_; // the final container (cellId -> faceId -> nodeId) for recosntruct a new mesh with FSDM
  const std::vector<FSFaceMatch>& matches_;
  const FSCell2NodeBuilder& cell2NodeBuilder_;
  FS_floatT tol_;



  FS_intT FindNodeLocalElem(const FSVec3& p,
                            const Cell2NodeData& cell2Node) const;
};

_FS_END_NAMESPACE
#endif
