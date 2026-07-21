#ifndef FSPOLYFACEBUILDER_H
#define FSPOLYFACEBUILDER_H

#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSCell2NodeBuilder.h"
#include <unordered_set>

_FS_BEGIN_NAMESPACE

struct FaceData {
  std::vector<FS_intT> nodeIds;
};


class FSPolyFaceBuilder
{
public:
  FSPolyFaceBuilder(const std::vector<FSFaceMatch>& matches, const FSCell2NodeBuilder& cell2NodeBuilder,
                    const std::unordered_set<GeomFaceKey, GeomFaceKeyHash>& faceKeys, FS_floatT tol = 1e-12)
    : matches_(matches), cell2NodeBuilder_(cell2NodeBuilder), boundaryFaceKeys_(faceKeys), tol_(tol) {};

  void AddInnerFaces(FSMeshEnums::CellType, const FSCellPool& cellPool, FS_intT cell, FS_intT face,
                     const FSFloatArrayT& oldCoords);

  void Build(FSMeshPolyFaceStorage& polyFaces);

  void CollectMatchesFaces();

  void Reorienting();

  FS_intT LocalCellIndex(FS_intT globalId) const;

  //! Whether a cell (global id) is known to the builder, i.e. it carries a
  //! clipped-surface match. Boundary cells without a match (possible in
  //! parallel: the local marker face has no counterpart in the other mesh)
  //! are absent and must not have inner faces built.
  bool HasCell(FS_intT globalId) const;

  const std::vector<std::vector<FaceData> >& CellFaces() const noexcept { return cellFaces_; };

  FSIntRegisterT Cell2NodePoly2D();

private:
  std::vector<std::vector<FaceData> > cellFaces_; // cellId → faceId → nodeId, used to reconstruct the new mesh
  const std::vector<FSFaceMatch>& matches_;
  const FSCell2NodeBuilder& cell2NodeBuilder_;
  const std::unordered_set<GeomFaceKey, GeomFaceKeyHash>& boundaryFaceKeys_;

  std::unordered_map<FS_intT, FS_intT> cellId2L_; // cellId → local index, for O(1) lookup
  FS_floatT tol_;

  FS_intT FindNodeLocalElem(const FSVec3& p, const Cell2NodeData& cell2Node) const;

  FSVec3 ComputeFaceCenter(const FaceData&, const std::vector<FSVec3>&) const;
  FSVec3 ComputeFaceNormal(const FaceData& face, const std::vector<FSVec3>& coords) const;
  FSVec3 ComputeCellCenter(const Cell2NodeData&) const;
};

_FS_END_NAMESPACE
#endif
