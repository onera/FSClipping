#include "FSClipping/FSTopologyBuilder.h"
#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSClipping/FSPolyFaceBuilder.h"

_FS_BEGIN_NAMESPACE

FSTopologyData FSTopologyBuilder::BuildSurfaceTopo(const std::vector<FSFaceMatch>& matches)
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

    cell2NodeBuilder_->AddCellNodes(elemIndex, cellType);

    if(f.type != FSFaceMatch::UNKNOWN) {
      cell2NodeBuilder_->AddClippedPolygon(elemIndex, f.clippedPoly3D);
    }
  }
  // cell2NodeBuilder_->UpdateNumCellsBorder();
  cell2NodeBuilder_->BuildGlobalNumbering();

  result.globalCoords = cell2NodeBuilder_->GlobalCoords();
  result.cell2Node = cell2NodeBuilder_->Cell2Node();
  // -------------------------------------------------
  // 2) Build polygonal faces
  // -------------------------------------------------
  // FSPolyFaceBuilder polyFaceBuilder(matches, *cell2NodeBuilder_, tol_);
  polyFaceBuilder_ = std::make_unique<FSPolyFaceBuilder>(FSPolyFaceBuilder(matches, *cell2NodeBuilder_, tol_));
  polyFaceBuilder_->CollectMatchesFaces();
  polyFaceBuilder_->Build(result.polyFaces);

  return result;
}


FSTopologyData FSTopologyBuilder::BuildVolumeTopo(const FSIntArrayT& cell2Node,
                                                  const std::set<FS_intT>& bdryCells,
                                                  const FSCellPool& cellPool,
                                                  const FSFloatArrayT& oldCoords)
{
  if(cell2NodeBuilder_->Cell2Node().Size() == 0)
    FSError.SetAndPrintAndExit("FSTopologyBuilder::BuildVolumeTopo : You need to build the surface topology before build the volume one.");

  FSTopologyData volumeResult;

  // FSConformalVolumeBuilder conformalVolumeBuilder(*cell2NodeBuilder_, tol_);
  // conformalVolumeBuilder.BuildVolume(cell2Node, bdryCells, type, oldCoord);

  auto type = cellPool.GetCellType();
  for(const auto& c : bdryCells)
    cell2NodeBuilder_->AddVolumeCellNodes(c, type, cell2Node, oldCoords);

  cell2NodeBuilder_->BuildGlobalNumbering();
  volumeResult.globalCoords = cell2NodeBuilder_->GlobalCoords();
  volumeResult.cell2Node = cell2NodeBuilder_->Cell2Node();

  // auto offSet = cellPool.GetOffset();
  // auto numCell = cellPool.GetNCells();
  // std::cout << "Offset du cell Pool : ";
  // for(FS_intT c = offSet; c < numCell + offSet; c++) {
  //   std::cout << c << " ";
  //   if(!bdryCells.contains(c))
  //     cell2NodeBuilder_->AddVolumeCellNodesInner(c, type, cell2Node, oldCoords); // pas besoin d'avoir deux fct AddVolume...Inner and AddVolumeCellNode
  // }
  // std::cout << "\n";



  for(const auto& c : bdryCells) {
    //        std::cout << "Element " << e;
    FS_intT nFaces = FSCellInfo::NFaces(type);
    for(FS_intT f = 0; f < nFaces; f++) {
      //        std::cout << "{ ";
      polyFaceBuilder_->AddInnerFaces(type, cellPool, c, f, oldCoords);

      //    std::cout << "} ";
    }
    // std::cout << "\n";
  }
  polyFaceBuilder_->Build(volumeResult.polyFaces);

  return volumeResult;
}
_FS_BEGIN_NAMESPACE