#include "FSClipping/FSTopologyBuilder.h"
#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSClipping/FSPolyFaceBuilder.h"

_FS_BEGIN_NAMESPACE

FSTopologyData FSTopologyBuilder::Build(const std::vector<FSFaceMatch>& matches)
{
  if(matches.empty())
    FSError.SetAndPrintAndExit("There was an error with the clipping algorithm. No matches were found. Look at the boundary marker for each mesh.");

  FSTopologyData result;

  // -------------------------------------------------
  // 1) Build cell2node + global coordinates
  // -------------------------------------------------
  FSCell2NodeBuilder cell2NodeBuilder(tol_);

  for(const auto& f : matches) {
    const auto elemIndex = f.elemOwner1;
    const auto cellType = f.elemOwnerType1;

    cell2NodeBuilder.AddOldCellNodes(elemIndex, cellType);

    if(f.type != FSFaceMatch::UNKNOWN) {
      cell2NodeBuilder.AddClippedPolygon(elemIndex, f.clippedPoly3D);
    }
  }

  cell2NodeBuilder.BuildGlobalNumbering();

  result.globalCoords = cell2NodeBuilder.GlobalCoords();
  result.cell2Node = cell2NodeBuilder.Cell2Node();

  // -------------------------------------------------
  // 2) Build polygonal faces
  // -------------------------------------------------
  FSPolyFaceBuilder polyFaceBuilder(matches, cell2NodeBuilder, tol_);
  polyFaceBuilder.Build(result.polyFaces);

  return result;
}

_FS_BEGIN_NAMESPACE