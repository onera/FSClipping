#ifndef FSNODEREGISTRY_H
#define FSNODEREGISTRY_H

#include "FSClipping/FSCellClipped.h"
#include "FSClipping/FSCellSource.h"
#include <FSVec3.h>
#include <span>
#include <vector>

_FS_BEGIN_NAMESPACE

using NodeId = FS_intT; // [0, nOld) : old node, same number ; [nOld, ...) : new polygon vertex

class FSNodeRegistry
{
public:
  FSNodeRegistry(FSMesh& mesh, const std::vector<FSCellSource>& sources, const FSCellClipped& clipped, FS_floatT tol);

  NodeId Old(FS_intT oldNode) const noexcept { return oldNode; } // identity, no tolerance
  std::span<const NodeId> Polygon(FS_intT match) const;
  const FSVec3& Coord(NodeId n) const;
  FS_intT ClipperId(NodeId n) const;
  // Exact identity shared by every proc: global number of an old node, nOldGlobal + ClipperId for a new vertex
  FS_intT GlobalId(NodeId n) const;
  FS_intT Size() const noexcept { return FS_intT(coords_.size()); }

private:
  std::vector<FSVec3> coords_;       // coordinates per NodeId
  std::vector<FS_intT> clipperId_;   // nodeGlobalIds per NodeId, -1 if the node is on no polygon
  std::vector<NodeId> polyNodes_;    // vertices of all the matches, flat
  std::vector<FS_intT> polyOffsets_; // match m -> [polyOffsets_[m], polyOffsets_[m+1]) in polyNodes_

  const FSCellPool* nodePool_ = nullptr; // old nodes, owned then ghost
  FS_intT nOld_ = 0;
  FS_intT nOldGlobal_ = 0;
};

_FS_END_NAMESPACE
#endif
