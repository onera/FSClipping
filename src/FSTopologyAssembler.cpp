#include "FSClipping/FSTopologyAssembler.h"
#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSClipping/FSPolyFaceBuilder.h"

_FS_BEGIN_NAMESPACE

FSTopologyData FSTopologyAssembler::BuildSurfaceTopo(const std::vector<FSFaceMatch>& matches)
{
  if(matches.empty())
    FSError.SetAndPrintAndExit("There was an error with the clipping algorithm. No matches were found. Look at the boundary marker for each mesh.");

  FSTopologyData result;

  // -------------------------------------------------
  // 1) Build cell2node + global coordinates
  // -------------------------------------------------
  // FSCell2NodeBuilder cell2NodeBuilder(tol_);

  for(const auto& f : matches) {
    const auto elemIndex = f.elemOwner1;
    const auto cellType = f.elemOwnerType1;

    cell2NodeBuilder_.AddCellNodes(elemIndex, cellType);

    if(f.type != FSFaceMatch::UNKNOWN) {
      cell2NodeBuilder_.AddClippedPolygon(elemIndex, f.clippedPoly3D);
    }
  }
  cell2NodeBuilder_.SetAllCellOnTheBorder();
  cell2NodeBuilder_.BuildGlobalNumbering();

  result.globalCoords = cell2NodeBuilder_.GlobalCoords();
  result.cell2NodePoly3D = cell2NodeBuilder_.Cell2NodePoly();
  // -------------------------------------------------
  // 2) Build polygonal faces
  // -------------------------------------------------
  polyFaceBuilder_ = std::make_unique<FSPolyFaceBuilder>(FSPolyFaceBuilder(matches, cell2NodeBuilder_, tol_));
  polyFaceBuilder_->CollectMatchesFaces();
  polyFaceBuilder_->Build(result.polyFaces);
  surfaceBuilt_ = true;

  return result;
}

void FSTopologyAssembler::CheckSurfaceWasBuilt() const
{
  if(!surfaceBuilt_)
    FSError.SetAndPrintAndExit("Surface topology must be built before volume topology.");
}

FSTopologyData FSTopologyAssembler::BuildVolumeTopo(const FSIntArrayT& cell2Node,
                                                    const std::set<FS_intT>& bdryCellPool,
                                                    const FSCellPool& cellPool,
                                                    const FSFloatArrayT& oldCoords)
{

  CheckSurfaceWasBuilt();

  FSTopologyData volumeResult;

  auto type = cellPool.GetCellType();

  auto offSet = cellPool.GetOffset();
  auto numCell = cellPool.GetNCells();
  for(FS_intT c = offSet; c < numCell + offSet; c++) {
    cell2NodeBuilder_.AddVolumeCellNodes(c, type, cell2Node, oldCoords);
  }
  cell2NodeBuilder_.BuildGlobalNumbering();

  volumeResult.globalCoords = cell2NodeBuilder_.GlobalCoords();
  volumeResult.cell2NodePoly3D = cell2NodeBuilder_.Cell2NodePoly();
  auto [it, inserted] = volumeResult.cell2NodeInner.try_emplace(type);
  it->second = cell2NodeBuilder_.Cell2NodeInner(type);

  for(const auto& c : bdryCellPool) {
    FS_intT nFaces = FSCellInfo::NFaces(type);
    for(FS_intT f = 0; f < nFaces; f++) {
      polyFaceBuilder_->AddInnerFaces(type, cellPool, c, f, oldCoords);
    }
  }
  polyFaceBuilder_->Build(volumeResult.polyFaces);


  return volumeResult;
}

// FSTopologyData FSTopologyAssembler::BuildVolumeTopo(const FSIntArrayT& cell2Node,
//                                                     const std::set<FS_intT>& bdryCellPool,
//                                                     const FSCellPool& cellPool,
//                                                     const FSFloatArrayT& oldCoords)
//{
//   CheckSurfaceWasBuilt();
//
//   FSTopologyData volumeResult;
//
//   auto type = cellPool.GetCellType();
//   for(const auto& c : bdryCellPool)
//     cell2NodeBuilder_.AddVolumeCellNodes(c, type, cell2Node, oldCoords);
//
//   cell2NodeBuilder_.BuildGlobalNumbering();
//   volumeResult.cell2NodePoly3D = cell2NodeBuilder_.Cell2NodePoly();
//
//   for(const auto& c : bdryCellPool) {
//     FS_intT nFaces = FSCellInfo::NFaces(type);
//     for(FS_intT f = 0; f < nFaces; f++) {
//       polyFaceBuilder_->AddInnerFaces(type, cellPool, c, f, oldCoords);
//     }
//   }
//   polyFaceBuilder_->Build(volumeResult.polyFaces);
//
//   auto offSet = cellPool.GetOffset();
//   auto numCell = cellPool.GetNCells();
//   for(FS_intT c = offSet; c < numCell + offSet; c++) {
//     if(!bdryCellPool.contains(c))
//       cell2NodeBuilder_.AddVolumeCellNodes(c, type, cell2Node, oldCoords);
//   }
//
//   cell2NodeBuilder_.BuildGlobalNumbering();
//   volumeResult.globalCoords = cell2NodeBuilder_.GlobalCoords();
//   auto [it, inserted] = volumeResult.cell2NodeInner.try_emplace(type);
//   it->second = cell2NodeBuilder_.Cell2NodeInner(type);
//
//   return volumeResult;
// }



_FS_BEGIN_NAMESPACE