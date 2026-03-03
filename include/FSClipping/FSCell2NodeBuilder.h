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
// We need that part because we have issue with std::unordered_map<NodeKey, ...>. If not, we have the following error :
// "Call to implicitly-deleted default constructor of 'std::unordered_map<NodeKey, FS_intT>'
// (aka 'unordered_map<NodeKey, int>')clang(ovl_deleted_special_init) FSCell2NodeBuilder.h(47, 40): Default constructed field 'coordToNode_' declared here"
// see : https://medium.com/@gulshansharma014/call-to-implicitly-deleted-default-constructor-of-unordered-map-pair-int-int-int-d3b2a6da0b41 or ChatGPT :P
// You can ignore the following part
namespace std {
template<>
struct hash<NodeKey> {
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
// and stop ignoring.

struct Cell2NodeData {
  std::vector<FSVec3> coords;                      // cell2Node with coordinates
  std::vector<FS_intT> nodeIds;                    // cell2node with local index in the globalCoords container
  std::unordered_map<NodeKey, FS_intT> localIndex; // mapping coord -> local index in coords vector. Very important later for the face2NodeBuilder

  bool onTheBorder = false;
};

/* --------------------------------------------------------------
   FSCell2NodeBuilder :
   This class reconstruct the connectivity cell2Node for all the
   cells.
   - First the methods AddVolumeCellNodes and AddClippedPolygon
   are used for add only the cell2Node.coords.
   - Second, in the BuildGlobalNumbering we store all the coords
   (without duplicate) in globalCoords_ and build the
   cell2Node.nodeIds (local index in globalCoords). Then we store
   cell2Node.localIndex, the mapping : coords -> local index in
   cell2Node (usefull later for face2Node).
   -------------------------------------------------------------- */
class FSCell2NodeBuilder
{
public:
  explicit FSCell2NodeBuilder(FS_floatT tol = 1e-12) : tol_(tol) {};

  void AddCellNodes(FS_intT cellId,
                    FSMeshEnums::CellType cellType);

  void AddVolumeCellNodesInner(FS_intT cellId,
                               FSMeshEnums::CellType cellType, // on peut se passer du cellType car cette classe est appele pour un type
                               const FSIntArrayT& cell2Node,
                               const FSFloatArrayT& coords);

  void AddVolumeCellNodes(FS_intT cellId,
                          FSMeshEnums::CellType cellType,
                          const FSIntArrayT& cell2Node,
                          const FSFloatArrayT& coords);

  void AddClippedPolygon(FS_intT cellId,
                         const std::vector<FSVec3>& poly);

  void BuildGlobalNumbering(); // reconstruct the global numbering for a new mesh

  void SetAllCellOnTheBorder();

  const std::unordered_map<FS_intT, Cell2NodeData>& CellData() const noexcept { return cellData_; };

  const std::vector<FSVec3>& GlobalCoords() const noexcept { return globalCoords_; };

  const std::unordered_map<NodeKey, FS_intT>& CoordToNode() const noexcept { return coordToNode_; };

  const std::unordered_map<FS_intT, FS_intT> CellId2L() const noexcept { return cellId2L_; };

  FSIntRegisterT Cell2NodePoly();

  FSIntArrayT Cell2NodeInner(FSMeshEnums::CellType type);

  FS_intT LocalCellIndex(FS_intT globalId) const;

  FS_intT GlobalCellId(FS_intT localId) const;

private:
  FS_floatT tol_;
  FS_intT numCellOnTheBorder_ = 0;
  std::unordered_map<FS_intT, Cell2NodeData> cellData_; // the main data with the mapping IdCellInFSDM -> coords, local node Id, local node Id
  std::vector<FS_intT> cellIds_;                        // sorted containers with the cell id
  std::unordered_map<FS_intT, FS_intT> cellId2L_;       // connectivity cellIds_[i] -> i, usefull for searching in O(1)

  /* Two following containers are global and use for searching if a node already exist in the border */
  std::vector<FSVec3> globalCoords_;                 // the list of all the coords of the border
  std::unordered_map<NodeKey, FS_intT> coordToNode_; // unicity mapping of each coord : coord -> local index in globalCood.

  bool exists(const std::vector<FSVec3>& nodes,
              const FSVec3& p) const;

  void BuildCellIdMapping(); // build the containers cellIds and cellId2L
};

_FS_END_NAMESPACE

#endif