#include "FSClipping/FSNodeRegistry.h"
#include <FSMeshData.h>
#include <cassert>
#include <unordered_map>

_FS_BEGIN_NAMESPACE

FSNodeRegistry::FSNodeRegistry(FSMesh& mesh, const std::vector<FSCellSource>& sources, const FSCellClipped& clipped,
                               FS_floatT tol)
{
  // 1 - Add the old coordinates
  FSFloatArrayT oldCoords;
  FS_intT nodeOffset = 0;
  FSQuantityDescArrayT coordDesc;

  if(!mesh.GetMeshData()->GetUnstructCells().GetCoordinates3D(coordDesc, oldCoords, nodeOffset))
    FSError.SetAndPrintAndExit("FSNodeRegistry: cannot get the node coordinates.");
  if(nodeOffset != 0)
    FSError.SetAndPrintAndExit("FSNodeRegistry: local node numbering expected (nodeOffset == 0).");

  nOld_ = oldCoords.Size(0);
  nodePool_ = mesh.GetMeshData()->GetUnstructCells().GetCellPool(FSMeshEnums::CT_Node);
  if(!nodePool_)
    FSError.SetAndPrintAndExit("FSNodeRegistry: the mesh has no node pool.");
  nOldGlobal_ = nodePool_->GetNGlobalCells();
  coords_.reserve(nOld_);
  for(FS_intT n = 0; n < nOld_; n++)
    coords_.emplace_back(oldCoords(n, 0), oldCoords(n, 1), oldCoords(n, 2));
  clipperId_.assign(nOld_, -1);

  // 2 - Groupe the cell source by type
  std::unordered_map<FS_intT, const FSCellSource*> sourceOfTypes;
  for(const auto& s : sources)
    sourceOfTypes.emplace(s.Type(), &s);

  // 3 - Complete the containers of the class
  std::unordered_map<FS_intT, NodeId> clipperToNode; // nodeGlobalIds -> NodeId, each vertex resolved once
  polyOffsets_.push_back(0);
  clipped.ForEachMatch([&](FS_intT, const FSFaceMatch& match) {
    const auto itSrc = sourceOfTypes.find(match.elemOwnerType1);
    if(itSrc == sourceOfTypes.end())
      FSError.SetAndPrintAndExit("FSNodeRegistry: no cell source for the owner type of a match.");

    const FSCellSource& src = *itSrc->second;
    const FS_intT cell = match.elemOwner1;
    const FS_intT nCellNodes = src.NNodes(cell);

    for(std::size_t i = 0; i < match.clippedPoly3D.size(); i++) {
      const FS_intT gId = match.nodeGlobalIds[i];
      const FSVec3& p = match.clippedPoly3D[i];

      // case 1: already resolved, exact clipper id, no tolerance
      const auto itNode = clipperToNode.find(gId);
      if(itNode != clipperToNode.end()) {
        polyNodes_.push_back(itNode->second);
        continue;
      }

      // case 2: coincides with a nodes in the old coord
      NodeId id = -1;
      FS_floatT best = tol;

      for(FS_intT n = 0; n < nCellNodes; n++) {
        const FS_intT oldNodes = src.Node(cell, n);
        const FS_floatT d = coords_[oldNodes].Distance(p);
        if(d <= best) {
          best = d;
          id = oldNodes;
        }
      }

      // case 3: new nodes to add
      if(id < 0) {
        id = Size();
        coords_.push_back(p);
        clipperId_.push_back(gId);
      } else {
        // an old node must not receive two different clipper ids
        assert(clipperId_[id] == -1 || clipperId_[id] == gId);
        clipperId_[id] = gId;
      }

      clipperToNode.emplace(gId, id);
      polyNodes_.push_back(id);
    }

    polyOffsets_.push_back(FS_intT(polyNodes_.size()));
  });
}

std::span<const NodeId> FSNodeRegistry::Polygon(FS_intT match) const
{
  FS_intT begin = polyOffsets_[match];
  FS_intT size = polyOffsets_[match + 1] - begin;
  std::span<const NodeId> spanNodes(polyNodes_);
  return spanNodes.subspan(begin, size);
}

const FSVec3& FSNodeRegistry::Coord(NodeId n) const
{
  if(n < 0 || n >= Size()) {
    FSString msg = "FSNodeRegistry: Coord at ";
    msg.Add(n);
    msg += " doesn't exist";
    FSError.SetAndPrintAndExit(msg.c_str());
  }
  return coords_[n];
}

FS_intT FSNodeRegistry::GlobalId(NodeId n) const
{
  if(n >= nOld_)
    return nOldGlobal_ + ClipperId(n);
  if(nodePool_->GetNumbering() == FSMeshEnums::CN_Local)
    return nodePool_->Local2GlobalNumber(n);
  return n; // global numbering with nodeOffset == 0: the local number is the global one
}

FS_intT FSNodeRegistry::ClipperId(NodeId n) const
{
  if(n < 0 || n >= Size()) {
    FSString msg = "FSNodeRegistry: clipperId at ";
    msg.Add(n);
    msg += " doesn't exist";
    FSError.SetAndPrintAndExit(msg.c_str());
  }
  return clipperId_[n];
}

_FS_END_NAMESPACE
