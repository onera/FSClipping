#ifndef FSCELL2NODEBUILDER_H
#define FSCELL2NODEBUILDER_H

#include <FSArray.h>
#include <FSRegister.h>
#include <FSTypes.h>
#include <FSVec3.h>
#include "FSMeshEnums.h"

_FS_BEGIN_NAMESPACE

struct NodeKey {
  FS_intT ix, iy, iz;

  NodeKey(const FSVec3& p, const FS_floatT tol);
  bool operator==(const NodeKey&) const = default;
};

// std::hash specialization required for NodeKey to be used as an
// unordered_map key. Without it, the default constructor of
// unordered_map<NodeKey, ...> is implicitly deleted.
namespace std {
template<> struct hash<NodeKey> {
  size_t operator()(const NodeKey& k) const noexcept
  {
    size_t h1 = std::hash<FS_intT>{}(k.ix);
    size_t h2 = std::hash<FS_intT>{}(k.iy);
    size_t h3 = std::hash<FS_intT>{}(k.iz);

    size_t seed = h1;
    seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};
} // namespace std

struct Cell2NodeData {
  std::vector<FSVec3> coords;                      // 3D coordinates of each node
  std::vector<FS_intT> nodeIds;                    // local index in globalCoords for each node
  std::unordered_map<NodeKey, FS_intT> localIndex; // coord → local index within this cell (used by FSPolyFaceBuilder)

  bool onTheBorder = false;
};

/* --------------------------------------------------------------
   FSCell2NodeBuilder
   Reconstructs the cell-to-node connectivity for all cells.

   Usage:
   1. Call AddCellNodes / AddVolumeCellNodes / AddClippedPolygon to
      populate coords per cell.
   2. Call BuildGlobalNumbering to assign unique indices, populate
      globalCoords_, coordToNode_, and cell2Node.localIndex.
   -------------------------------------------------------------- */
class FSCell2NodeBuilder
{
public:
  explicit FSCell2NodeBuilder(FS_floatT tol = 1e-12) : tol_(tol) {};

  void AddCellNodes(FS_intT cellId, FSMeshEnums::CellType cellType);

  void AddVolumeCellNodesInner(FS_intT cellId, FSMeshEnums::CellType cellType, const FSIntArrayT& cell2Node,
                               const FSFloatArrayT& coords);

  void AddVolumeCellNodes(FS_intT cellId, FSMeshEnums::CellType cellType, const FSIntArrayT& cell2Node,
                          const FSFloatArrayT& coords);

  void AddClippedPolygon(FS_intT cellId, const std::vector<FSVec3>& poly);

  // Same as above but also records the clipper-assigned global node IDs
  // (level 2 parallelism). globalIds is parallel to poly.
  void AddClippedPolygon(FS_intT cellId, const std::vector<FSVec3>& poly, const std::vector<FS_intT>& globalIds);

  void BuildGlobalNumbering();

  // For each local node index (in GlobalCoords order), the clipper-assigned
  // global ID, or -1 if the node was never fed through the globalIds overload.
  // Requires BuildGlobalNumbering to have run.
  FSIntArrayT NodeGlobalNumbers() const;

  void SetAllCellOnTheBorder();

  const std::unordered_map<FS_intT, Cell2NodeData>& CellData() const noexcept { return cellData_; };

  const std::vector<FSVec3>& GlobalCoords() const noexcept { return globalCoords_; };

  const std::unordered_map<NodeKey, FS_intT>& CoordToNode() const noexcept { return coordToNode_; };

  const std::unordered_map<FS_intT, FS_intT> CellId2L() const noexcept { return cellId2L_; };

  const FS_intT NumCellOnTheBorder() const noexcept { return numCellOnTheBorder_; };

  FSIntRegisterT Cell2NodePoly3D();

  FSIntArrayT Cell2NodeInner(FSMeshEnums::CellType type);

  FS_intT LocalCellIndex(FS_intT globalId) const;

  FS_intT GlobalCellId(FS_intT localId) const;

private:
  FS_floatT tol_;
  FS_intT numCellOnTheBorder_ = 0;
  std::unordered_map<FS_intT, Cell2NodeData> cellData_; // FSDM cell id → node data
  std::vector<FS_intT> cellIds_;                        // sorted cell ids (parallel to local indices)
  std::unordered_map<FS_intT, FS_intT> cellId2L_;       // cellIds_[i] → i, for O(1) lookup

  std::vector<FSVec3> globalCoords_;                 // deduplicated list of all node coordinates
  std::unordered_map<NodeKey, FS_intT> coordToNode_; // coord → index in globalCoords_
  std::unordered_map<NodeKey, FS_intT> keyToGlobalId_; // coord → clipper-assigned global node ID

  bool exists(const std::vector<FSVec3>& nodes, const FSVec3& p) const;

  void BuildCellIdMapping();
};

_FS_END_NAMESPACE

#endif
