#include "FSClipping/FSTopologyAssembler.h"
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSClipping/FSPolyFaceBuilder.h"
#include <FSMeshData.h>
#include <fstream>
#include <sstream>
#include <unistd.h>

_FS_BEGIN_NAMESPACE

FSTopologyData FSTopologyAssembler::BuildSurfaceTopo(const std::vector<FSFaceMatch>& matches,
                                                     const std::unordered_set<GeomFaceKey, GeomFaceKeyHash>& faceKeys)
{
  FSTopologyData result;

  // In parallel a proc whose sub-domain does not cut the interface has no
  // match. That is legitimate: build an empty (but valid) surface topology so
  // the proc still participates in the collective reconstruction. The volume
  // step then passes all its cells through as unclipped volume cells.
  if(matches.empty()) {
    polyFaceBuilder_ = std::make_unique<FSPolyFaceBuilder>(matches, cell2NodeBuilder_, faceKeys, tol_);
    polyFaceBuilder_->CollectMatchesFaces();
    surfaceBuilt_ = true;
    return result;
  }

  // -------------------------------------------------
  // DEBUG: Write elemOwner1 and elemOwner2 to per-process files
  // -------------------------------------------------
  // Uses getpid() for unique per-process file names in parallel runs.
  {
    pid_t pid = getpid();
    std::ostringstream oss;
    oss << "elemOwner_proc" << pid << ".txt";
    std::ofstream outFile(oss.str());
    if(outFile.is_open()) {
      for(const auto& f : matches) {
        outFile << "elemOwner1: " << f.elemOwner1 << " | elemOwner2: " << f.elemOwner2 << "\n";
        for(const auto& ng : f.nodeGlobalIds) {
          outFile << ng << ", ";
        }
        outFile << "\n";
      }
      outFile.close();
    }
  }
  // -------------------------------------------------
  // 1. Build cell2node + global coordinates
  // -------------------------------------------------
  for(const auto& f : matches) {
    const auto elemIndex = f.elemOwner1;
    const auto cellType = f.elemOwnerType1;

    cell2NodeBuilder_.AddCellNodes(elemIndex, cellType);

    if(f.type != FSFaceMatch::UNKNOWN) {
      if(!f.nodeGlobalIds.empty())
        cell2NodeBuilder_.AddClippedPolygon(elemIndex, f.clippedPoly3D, f.nodeGlobalIds);
      else
        cell2NodeBuilder_.AddClippedPolygon(elemIndex, f.clippedPoly3D);
    }
  }
  cell2NodeBuilder_.SetAllCellOnTheBorder();
  cell2NodeBuilder_.BuildGlobalNumbering();

  // -------------------------------------------------
  // DEBUG: Write Cell2NodeBuilder to per-process files
  // -------------------------------------------------
  // Uses getpid() for unique per-process file names in parallel runs.
  {
    pid_t pid = getpid();
    std::ostringstream oss;
    oss << "cellData_proc" << pid << ".txt";
    std::ofstream outFile(oss.str());
    if(outFile.is_open()) {
      for(const auto& elem : cell2NodeBuilder_.CellData()) {
        outFile << elem.first << " : ";
        for(const auto& c : elem.second.nodeIds) {
          NodeKey key(c, tol_);
          outFile << c << " ";
        }
        outFile << "\n";
      }
      outFile.close();
    }
  }
  result.globalCoords = cell2NodeBuilder_.GlobalCoords();
  result.cell2NodePoly3D = cell2NodeBuilder_.Cell2NodePoly3D();

  // Global numbering assigned by the clipper (level 2): expose it to
  // FSMeshReconstruction so that procs of the same mesh sub-communicator
  // agree on shared interface nodes. Absent in sequential mode.
  const bool hasGlobalIds = matches.front().globalCellId >= 0;
  if(hasGlobalIds) {
    result.nodeGlobalNumbers = cell2NodeBuilder_.NodeGlobalNumbers();
    for(const auto& f : matches)
      result.poly2DGlobalNumbers.Append(f.globalCellId);
  }
  // -------------------------------------------------
  // 2. Build polygonal faces
  // -------------------------------------------------
  polyFaceBuilder_ = std::make_unique<FSPolyFaceBuilder>(FSPolyFaceBuilder(matches, cell2NodeBuilder_, faceKeys, tol_));
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

void FSTopologyAssembler::BuildVolumeTopo(const FSIntArrayT& cell2Node, const std::set<FS_intT>& bdryCellPool,
                                          const FSCellPool& cellPool, const FSFloatArrayT& oldCoords,
                                          FSTopologyData& volumeResult)
{
  CheckSurfaceWasBuilt();

  auto type = cellPool.GetCellType();
  auto offSet = cellPool.GetOffset();
  auto nCells = cellPool.GetNCells();

  // -------------------------------------------------
  // 1. Add volume cells
  // -------------------------------------------------
  // Skip ghost cells: with several procs per mesh (local numbering), non-owned
  // cells are replicated on neighbour procs and would be emitted twice.
  for(FS_intT c = offSet; c < nCells + offSet; c++) {
    if(!cellPool.IsOwned(c))
      continue;
    cell2NodeBuilder_.AddVolumeCellNodes(c, type, cell2Node, oldCoords);
  }

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

  // Iterate the cells actually held by the builder (owned cells only — ghosts
  // were skipped above, so the pool count may exceed the builder count).
  const FS_intT nBuilderCells = static_cast<FS_intT>(cell2NodeBuilder_.CellData().size());
  for(FS_intT localIdx = 0; localIdx < nBuilderCells; ++localIdx) {
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
    // NOTE: if multiple element types touch the border, the last one wins here.
    volumeResult.cellParent[FSMeshEnums::CellType::CT_Poly3D] = std::move(cellParentPoly3D);
    volumeResult.cellParentType[FSMeshEnums::CellType::CT_Poly3D] = std::move(cellParentPolyType);
  }

  // -------------------------------------------------
  // 3. Build internal faces for clipped cells
  // -------------------------------------------------
  // Inner faces only for owned, matched boundary cells:
  //  - a ghost boundary cell is owned by another proc (skipped above), absent
  //    from the builder;
  //  - in parallel a cell may have a marker face locally but no clipping match
  //    (no counterpart in the other mesh), so it stays a plain volume cell with
  //    no Poly3D inner faces.
  // Either way it is absent from the poly-face builder → skip it.
  for(const auto& c : bdryCellPool) {
    if(!cellPool.IsOwned(c) || !polyFaceBuilder_->HasCell(c))
      continue;

    FS_intT nFaces = FSCellInfo::NFaces(type);
    for(FS_intT f = 0; f < nFaces; f++)
      polyFaceBuilder_->AddInnerFaces(type, cellPool, c, f, oldCoords);
  }
  polyFaceBuilder_->Reorienting();
  polyFaceBuilder_->Build(volumeResult.polyFaces);
  volumeResult.cell2NodePoly2D = polyFaceBuilder_->Cell2NodePoly2D();
}

FSIntArrayT FSTopologyAssembler::UpdateOldCell2Node(const FSIntArrayT& oldCell2Node, const FSCellPool& cellPool,
                                                    const FSFloatArrayT& oldCoords,
                                                    const std::set<FS_intT>& bdry2DCells, FSIntArrayT& cellParent)
{
  const FS_intT nCellOld = oldCell2Node.Size(0);
  const FS_intT nCellNodes = oldCell2Node.Size(1);
  const FS_intT offset = oldCell2Node.Offset();

  // Keep only owned, non-clipped surface cells. Ghost cells are replicated on
  // neighbour procs and their nodes were never registered in the builder
  // (BuildVolumeTopo skips ghosts) — emitting them would duplicate faces and
  // reference unknown nodes.
  FS_intT nKept = 0;
  for(FS_intT c = 0; c < nCellOld; c++) {
    const FS_intT cellId = offset + c;
    if(!bdry2DCells.contains(cellId) && cellPool.IsOwned(cellId))
      ++nKept;
  }

  FSIntArrayT cell2Node(nKept, nCellNodes);
  FS_intT iter = 0;

  for(FS_intT c = 0; c < nCellOld; c++) {
    const FS_intT cellId = offset + c;
    if(bdry2DCells.contains(cellId) || !cellPool.IsOwned(cellId))
      continue;

    cellParent.Append(cellId);
    for(FS_intT node = 0; node < nCellNodes; ++node) {
      FS_intT idx = oldCell2Node(cellId, node);
      FSVec3 vecNode{oldCoords(idx, 0), oldCoords(idx, 1), oldCoords(idx, 2)};
      // Register the node if it is unknown: an owned unclipped surface cell may
      // reference a node that no owned volume cell owns (shared only with volume
      // ghosts, which BuildVolumeTopo skips). topo.globalCoords is refreshed by
      // the caller after all surface types are processed.
      cell2Node(iter, node) = cell2NodeBuilder_.ResolveOrRegisterNode(vecNode);
    }
    iter++;
  }
  return cell2Node;
}

void FSTopologyAssembler::AppendUnclippedSurfaces(FSMesh& mesh,
                                                  const std::unordered_map<FS_intT, std::set<FS_intT> >& bdry2DCells,
                                                  const FSFloatArrayT& oldCoords, FSTopologyData& topo)
{
  for(const auto& t : mesh.GetCellTypes()) {
    if(!FSMeshEnums::IsUnstructSurfaceCellType(t))
      continue;

    const auto& oldCell2Node = mesh.GetCell2Node(t);
    const auto& cellPool = *mesh.GetMeshData()->GetUnstructCells().GetCellPool(t);

    const auto& bdry = bdry2DCells.contains(t) ? bdry2DCells.at(t) : std::set<FS_intT>{};

    FSIntArrayT parent;
    auto cell2Node = UpdateOldCell2Node(oldCell2Node, cellPool, oldCoords, bdry, parent);

    topo.cell2NodeInner[t] = std::move(cell2Node);
    topo.cellParent[t] = std::move(parent);
  }

  // UpdateOldCell2Node may have appended nodes that no owned volume cell owns.
  // Refresh globalCoords so the reconstructed mesh has coordinates for them.
  topo.globalCoords = cell2NodeBuilder_.GlobalCoords();
}

_FS_END_NAMESPACE
