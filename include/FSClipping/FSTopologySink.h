#ifndef FSTOPOLOGYSINK_H
#define FSTOPOLOGYSINK_H

#include "FSClipping/FSCellSource.h"
#include "FSClipping/FSNodeRegistry.h"
#include "FSClipping/FSPolyCell.h"
#include "FSClipping/FSTopologyData.h"
#include <algorithm>
#include <cassert>
#include <map>
#include <span>
#include <vector>

_FS_BEGIN_NAMESPACE

// NodeId -> index in the compacted output, -1 if the node is unused
struct FSRenumber {
  std::vector<FS_intT> newId;
  FS_intT operator()(NodeId n) const
  {
    assert(newId[n] >= 0);
    return newId[n];
  }
};

// One accumulator per FSDM container, filled during the loop and built once in Finalize
struct FSFixedRows {
  FS_intT nPerRow = 0;
  std::vector<FS_intT> values; // nPerRow NodeIds per row, flat

  FSIntArrayT Build(const FSRenumber&) const; // (nRows, nPerRow), values renumbered
};
struct FSRaggedRows {
  std::vector<FS_intT> counts, values; // counts[i] NodeIds for row i, flat
  FSIntRegisterT Build(const FSRenumber&) const; // values renumbered
};
struct FSPolyFaceRows {
  std::vector<FS_intT> nFaces, nFaceNodes, nodes; // nodes are positions in the cell row
  void Build(FSMeshPolyFaceStorage&) const; // local positions: never renumbered
};
struct FSParentRows {
  std::vector<FS_intT> parent, parentType; // parent cell (with pool offset) and its type
  void Build(FSIntArrayT& parentOut, FSIntArrayT& parentTypeOut) const;
};

// Keeps the used nodes in registry order, fills coords and clipperIds, returns NodeId -> new index
FSRenumber CompactNodes(const std::vector<bool>& used, const FSNodeRegistry& nodes, std::vector<FSVec3>& coords,
                        FSIntArrayT& clipperIds);

// Receives each output cell once, in NodeId; Finalize compacts the nodes and writes FSTopologyData
class FSTopologySink
{
public:
  explicit FSTopologySink(const FSNodeRegistry& nodes) : nodes_(nodes), used_(nodes.Size(), false) {}

  void AddStandard(const FSCellSource& src, FS_intT cell); // inner_[type] + parents_[type]
  void AddPoly2D(FS_intT parent, FSMeshEnums::CellType parentType, std::span<const NodeId> row,
                 FS_intT marker); // poly2D_ + parents_[Poly2D] + poly2DMarkers_
  void AddPoly3D(FS_intT parent, FSMeshEnums::CellType parentType,
                 const FSPolyCell&);                     // poly3D_ + faces_ + parents_[Poly3D]
  void AddPoly3D(const FSCellSource& src, FS_intT cell); // original Poly3D: row and faces copied

  FSTopologyData Finalize() &&; // call once: std::move(sink).Finalize()

private:
  // Row = nodes by first appearance in the faces; faces stored as positions in that row
  template<class NFaceNodes, class FaceNode>
  void AppendPolyhedron(FS_intT nFaces, NFaceNodes&& nFaceNodes, FaceNode&& faceNode)
  {
    assert(nFaces >= 4);
    faces_.nFaces.push_back(nFaces);
    for(FS_intT f = 0; f < nFaces; ++f) {
      const FS_intT nk = nFaceNodes(f);
      assert(nk >= 3);
      faces_.nFaceNodes.push_back(nk);
      for(FS_intT k = 0; k < nk; ++k) {
        const NodeId node = faceNode(f, k);
        // position of the node in the row, appended at its first appearance
        auto it = std::find(row_.begin(), row_.end(), node);
        if(it == row_.end()) {
          row_.push_back(node);
          it = row_.end() - 1;
        }
        faces_.nodes.push_back(static_cast<FS_intT>(it - row_.begin()));
      }
    }
    for(NodeId n : row_)
      used_[n] = true;
    poly3D_.counts.push_back(static_cast<FS_intT>(row_.size()));
    poly3D_.values.insert(poly3D_.values.end(), row_.begin(), row_.end());
    row_.clear();
  }

  const FSNodeRegistry& nodes_;
  std::vector<bool> used_;                                // globalCoords, nodeGlobalNumbers, FSRenumber
  std::map<FSMeshEnums::CellType, FSFixedRows> inner_;    // cell2NodeInner
  FSRaggedRows poly2D_, poly3D_;                          // cell2NodePoly2D, cell2NodePoly3D
  FSPolyFaceRows faces_;                                  // polyFaces
  std::map<FSMeshEnums::CellType, FSParentRows> parents_; // cellParent, cellParentType
  std::vector<FS_intT> poly2DMarkers_;                    // poly2DMarkers, aligned on poly2D_ (-1 = not carried)
  std::vector<NodeId> row_;                               // row of the cell being appended
};

_FS_END_NAMESPACE

#endif
