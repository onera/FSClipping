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
  result.cell2NodePoly3D = cell2NodeBuilder_.Cell2NodePoly3D();

  // -------------------------------------------------
  // 2) Build polygonal faces
  // -------------------------------------------------
  polyFaceBuilder_ = std::make_unique<FSPolyFaceBuilder>(FSPolyFaceBuilder(matches, cell2NodeBuilder_, tol_));
  polyFaceBuilder_->CollectMatchesFaces();
  polyFaceBuilder_->Build(result.polyFaces);

  result.cell2NodePoly2D = polyFaceBuilder_->Cell2NodePoly2D();
  FSIntArrayT cellParentPoly2D;
  FSIntArrayT cellParentPolyType;
  for(const auto& f : matches) {
    cellParentPoly2D.Append(f.faceOwner1);
    cellParentPolyType.Append(f.faceOwnerType1);
  }
  result.cellParent[FSMeshEnums::CellType::CT_Poly2D] = std::move(cellParentPoly2D);
  result.cellParentType[FSMeshEnums::CellType::CT_Poly2D] = std::move(cellParentPolyType);

  surfaceBuilt_ = true;

  return result;
}

void FSTopologyAssembler::CheckSurfaceWasBuilt() const
{
  if(!surfaceBuilt_)
    FSError.SetAndPrintAndExit("Surface topology must be built before volume topology.");
}

void FSTopologyAssembler::BuildVolumeTopo(const FSIntArrayT& cell2Node,
                                          const std::set<FS_intT>& bdryCellPool,
                                          const FSCellPool& cellPool,
                                          const FSFloatArrayT& oldCoords,
                                          FSTopologyData& volumeResult)
{
  CheckSurfaceWasBuilt();

  auto type = cellPool.GetCellType();

  auto offSet = cellPool.GetOffset();
  auto nCells = cellPool.GetNCells();

  // -------------------------------------------------
  // 1. Add volume cells
  // -------------------------------------------------

  for(FS_intT c = offSet; c < nCells + offSet; c++)
    cell2NodeBuilder_.AddVolumeCellNodes(c, type, cell2Node, oldCoords);

  cell2NodeBuilder_.BuildGlobalNumbering();

  // -------------------------------------------------
  // 2. Update topology data
  // -------------------------------------------------

  volumeResult.globalCoords = cell2NodeBuilder_.GlobalCoords();
  volumeResult.cell2NodePoly3D = cell2NodeBuilder_.Cell2NodePoly3D();

  auto& inner = volumeResult.cell2NodeInner[type];
  inner = cell2NodeBuilder_.Cell2NodeInner(type);

  FSIntArrayT cellParent;
  FSIntArrayT cellParentPolyType;
  FSIntArrayT cellParentPoly3D;

  for(FS_intT localIdx = 0; localIdx < nCells; ++localIdx) {
    FS_intT globalId = cell2NodeBuilder_.GlobalCellId(localIdx);
    const auto& data = cell2NodeBuilder_.CellData().at(globalId);
    if(data.onTheBorder) {
      cellParentPoly3D.Append(globalId);
      cellParentPolyType.Append(type);
    } else
      cellParent.Append(globalId);
  }

  volumeResult.cellParent[type] = std::move(cellParent);
  if(!cellParentPoly3D.IsEmpty()) {
    volumeResult.cellParent[FSMeshEnums::CellType::CT_Poly3D] = std::move(cellParentPoly3D); // attention si on a different type d element sur la frontiere, on ecrase
    volumeResult.cellParentType[FSMeshEnums::CellType::CT_Poly3D] = std::move(cellParentPolyType);
  }
  // -------------------------------------------------
  // 3. Build internal faces for clipped cells
  // -------------------------------------------------

  for(const auto& c : bdryCellPool) {
    FS_intT nFaces = FSCellInfo::NFaces(type);

    for(FS_intT f = 0; f < nFaces; f++)
      polyFaceBuilder_->AddInnerFaces(type, cellPool, c, f, oldCoords);
  }

  polyFaceBuilder_->Build(volumeResult.polyFaces);
  volumeResult.cell2NodePoly2D = polyFaceBuilder_->Cell2NodePoly2D();
}


FSIntArrayT FSTopologyAssembler::UpdateOldCell2Node(const FSIntArrayT& oldCell2Node,
                                                    const FSFloatArrayT& oldCoords,
                                                    const std::set<FS_intT>& bdry2DCells,
                                                    FSIntArrayT& cellParent)
{
  FS_intT nCellOld = oldCell2Node.Size(0);
  FS_intT nCellNew = nCellOld - bdry2DCells.size();
  FS_intT nCellNodes = oldCell2Node.Size(1);

  FSIntArrayT cell2Node(nCellNew, nCellNodes);
  const auto& coord2Node = cell2NodeBuilder_.CoordToNode();
  FS_intT offset = oldCell2Node.Offset();
  FS_intT iter = 0;

  for(FS_intT c = 0; c < nCellOld; c++) {
    if(!bdry2DCells.contains(offset)) {
      cellParent.Append(offset);
      for(FS_intT node = 0; node < nCellNodes; ++node) {
        FS_intT idx = oldCell2Node(offset, node);
        auto x = oldCoords(idx, 0);
        auto y = oldCoords(idx, 1);
        auto z = oldCoords(idx, 2);
        FSVec3 vecNode{x, y, z};
        NodeKey key(vecNode, tol_);
        cell2Node(iter, node) = coord2Node.at(key);
      }
      iter++;
    }
    offset++;
  }
  return cell2Node;
}

void FSTopologyAssembler::AppendUnclippedSurfaces(
  FSMesh& mesh,
  const std::unordered_map<FS_intT, std::set<FS_intT> >& bdry2DCells,
  const FSFloatArrayT& oldCoords,
  FSTopologyData& topo)
{
  for(const auto& t : mesh.GetCellTypes()) {
    if(!FSMeshEnums::IsUnstructSurfaceCellType(t))
      continue;

    const auto& oldCell2Node = mesh.GetCell2Node(t);

    const auto& bdry =
      bdry2DCells.contains(t) ? bdry2DCells.at(t)
                              : std::set<FS_intT>{};

    FSIntArrayT parent;

    auto cell2Node = UpdateOldCell2Node(oldCell2Node,
                                        oldCoords,
                                        bdry,
                                        parent);

    topo.cell2NodeInner[t] = std::move(cell2Node);
    topo.cellParent[t] = std::move(parent);
  }
}

_FS_BEGIN_NAMESPACE