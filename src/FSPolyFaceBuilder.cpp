#include "FSClipping/FSPolyFaceBuilder.h"

_FS_BEGIN_NAMESPACE

void FSPolyFaceBuilder::CollectMatchesFaces()
{
  const auto& cellData = cell2NodeBuilder_.CellData();
  cellFaces_.resize(cellData.size());

  for(const auto& f : matches_) {
    // std::cout << f.type << std::endl;
    if(f.type == FSFaceMatch::UNKNOWN)
      continue;

    FS_intT cellId = cell2NodeBuilder_.LocalCellIndex(f.elemOwner1);

    FaceData face;
    // std::cout << "Match : " << f.face1 << " with " << f.face2 << " element : " << f.elemOwner1 << "\n";
    // std::cout << "Type de match : " << f.type << std::endl;
    for(const auto& p : f.clippedPoly3D) {
      // auto tab = cellData.at(f.elemOwner1).coords;
      // std::cout << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")" << " ";
      face.nodeIds.push_back(FindNodeLocalElem(p, cellData.at(f.elemOwner1)));
      // std::cout << "\n";
    }
    // std::cout << "\n"
    //           << "\n";

    cellFaces_[cellId].push_back(std::move(face));
  }
}

FS_intT FSPolyFaceBuilder::FindNodeLocalElem(const FSVec3& p,
                                             const Cell2NodeData& cellData) const
{
  NodeKey k(p, tol_);

  auto it = cellData.localIndex.find(k); // searching in O(1), that's the reason of cell2Node.localIndex mapping
  if(it == cellData.localIndex.end()) {
    FSString msg = "Node : (";
    msg += std::to_string(p[0]);
    msg += ", ";
    msg += std::to_string(p[1]);
    msg += ", ";
    msg += std::to_string(p[2]);
    msg += ")";
    msg += " not found in the face : \n";
    for(const auto& v : cellData.coords) {
      msg += "(" + std::to_string(v[0]) + ", " +
             std::to_string(v[1]) + ", " +
             std::to_string(v[2]) + ") ";
    }
    msg += "\n";
    FSError.SetAndPrintAndExit(msg.c_str());
  }

  return it->second;
}

void FSPolyFaceBuilder::AddInnerFaces(FSMeshEnums::CellType type,
                                      const FSCellPool& cellPool,
                                      FS_intT cell,
                                      FS_intT face,
                                      const FSFloatArrayT& oldCoords)
{
  if(cellFaces_.empty())
    FSError.SetAndPrintAndExit("FSPolyFaceBuilder::AddInnerFaces cellFaces_ is empty. It's must be build during the surface topo build");

  FS_intT cellId = cell2NodeBuilder_.LocalCellIndex(cell);
  const auto& cellData = cell2NodeBuilder_.CellData();
  // FS_intT nFacesSurfaces = cellFaces_[localCellId].size();
  FS_intT nCorners = FSCellInfo::NFaceCorners(type, cellPool, cell, face);
  FaceData innerFace;
  innerFace.nodeIds.resize(nCorners);
  for(FS_intT fc = 0; fc < nCorners; fc++) {
    FS_intT node = FSCellInfo::GetCellFaceCorner(type, cellPool, cell, face, fc);
    FSVec3 p(oldCoords(node, 0), oldCoords(node, 1), oldCoords(node, 2));
    innerFace.nodeIds[fc] = FindNodeLocalElem(p, cellData.at(cell));
    // face.nodeIds.push_back(FindNodeLocalElem(p, cellData.at(f.elemOwner1)));
  }
  cellFaces_[cellId].push_back(std::move(innerFace));
}

void FSPolyFaceBuilder::Build(FSMeshPolyFaceStorage& polyFaces)
{
  // CollectMatchesFaces();

  FS_intT nCells = cellFaces_.size();
  polyFaces.Init(nCells);

  // 1-Collect the number of faces per poly
  for(FS_intT c = 0; c < nCells; ++c) {
    polyFaces.SetNumFacesForPoly(c, cellFaces_[c].size());
  }
  polyFaces.EndNumFacesSetup();

  // 2-Collect the number of nodes per faces
  for(FS_intT c = 0; c < nCells; ++c) {
    // std::cout << "Element " << c << " { ";
    for(std::size_t f = 0; f < cellFaces_[c].size(); ++f) {
      // std::cout << cellFaces_[c][f].nodeIds.size() << " ";
      polyFaces.SetNumNodesForFace(
        c, f, cellFaces_[c][f].nodeIds.size());
    }
    // std::cout << "}" << std::endl;
  }
  polyFaces.EndNumNodesSetup();

  // 3-Collect the local index of each coord in cell2Node.coords
  for(FS_intT c = 0; c < nCells; ++c) {
    // std::cout << "Element " << c << " : ";
    for(std::size_t f = 0; f < cellFaces_[c].size(); ++f) {
      // std::cout << " { ";
      for(FS_intT nid : cellFaces_[c][f].nodeIds) {
        polyFaces.AddFaceNode(c, f, nid);
        // std::cout << nid << " ";
      }
      // std::cout << "}";
    }
    // std::cout << "\n";
  }
}
_FS_END_NAMESPACE
