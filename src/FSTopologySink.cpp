#include "FSClipping/FSTopologySink.h"

_FS_BEGIN_NAMESPACE

void FSTopologySink::AddStandard(const FSCellSource& src, FS_intT cell)
{
  // 1 - Add the row to inner_[type]
  const auto type = src.Type();
  const FS_intT nNodes = src.NNodes(cell);

  auto [it, inserted] = inner_.try_emplace(type);
  FSFixedRows& rows = it->second;
  if(inserted)
    rows.nPerRow = nNodes;
  assert(rows.nPerRow == nNodes);

  for(FS_intT k = 0; k < nNodes; ++k) {
    const NodeId n = nodes_.Old(src.Node(cell, k));
    rows.values.push_back(n);
    used_[n] = true;
  }

  // 2 - Complete parents_[type]
  auto& p = parents_[type];
  p.parent.push_back(cell);
  p.parentType.push_back(type);
}

void FSTopologySink::AddPoly3D(const FSCellSource& src, FS_intT cell)
{
  // Row seeded with the original one, so the local face indices are kept
  row_ = src.Row(cell);
  [[maybe_unused]] const std::size_t nRow = row_.size();

  AppendPolyhedron(
    src.NFaces(cell), [&](FS_intT f) { return src.NFaceCorners(cell, f); },
    [&](FS_intT f, FS_intT k) { return src.FaceCorner(cell, f, k); });
  assert(poly3D_.counts.back() == static_cast<FS_intT>(nRow));

  auto& p = parents_[FSMeshEnums::CT_Poly3D];
  p.parent.push_back(cell);
  p.parentType.push_back(FSMeshEnums::CT_Poly3D);
}

void FSTopologySink::AddPoly3D(FS_intT parent, FSMeshEnums::CellType parentType, const FSPolyCell& poly)
{
  assert(row_.empty());
  AppendPolyhedron(
    static_cast<FS_intT>(poly.faces.size()), [&](FS_intT f) { return static_cast<FS_intT>(poly.faces[f].size()); },
    [&](FS_intT f, FS_intT k) { return poly.faces[f][k]; });

  auto& p = parents_[FSMeshEnums::CT_Poly3D];
  p.parent.push_back(parent);
  p.parentType.push_back(parentType);
}

void FSTopologySink::AddPoly2D(FS_intT parent, FSMeshEnums::CellType parentType, std::span<const NodeId> row,
                               FS_intT marker)
{
  assert(row.size() >= 3);
  poly2D_.counts.push_back(static_cast<FS_intT>(row.size()));
  for(NodeId n : row) {
    poly2D_.values.push_back(n);
    used_[n] = true;
  }

  auto& p = parents_[FSMeshEnums::CT_Poly2D];
  p.parent.push_back(parent);
  p.parentType.push_back(parentType);
  poly2DMarkers_.push_back(marker);
}

FSRenumber CompactNodes(const std::vector<bool>& used, const FSNodeRegistry& nodes, std::vector<FSVec3>& coords,
                        FSIntArrayT& clipperIds)
{
  assert(FS_intT(used.size()) == nodes.Size());
  const FS_intT nUsed = FS_intT(std::count(used.begin(), used.end(), true));

  FSRenumber renumber{std::vector<FS_intT>(used.size(), -1)};
  coords.clear();
  coords.reserve(nUsed);
  clipperIds = FSIntArrayT(nUsed);

  for(NodeId n = 0; n < nodes.Size(); ++n) {
    if(!used[n])
      continue;
    const FS_intT k = FS_intT(coords.size());
    renumber.newId[n] = k;
    coords.push_back(nodes.Coord(n));
    clipperIds[k] = nodes.ClipperId(n);
  }
  return renumber;
}

namespace {

FSIntArrayT ToIntArray(const std::vector<FS_intT>& v)
{
  FSIntArrayT a;
  if(v.empty())
    return a;
  a.Resize(FS_intT(v.size()));
  for(FS_intT i = 0; i < FS_intT(v.size()); ++i)
    a[i] = v[i];
  return a;
}

} // namespace

FSIntArrayT FSFixedRows::Build(const FSRenumber& renumber) const
{
  FSIntArrayT rows;
  if(values.empty())
    return rows;
  assert(nPerRow > 0 && FS_intT(values.size()) % nPerRow == 0);

  const FS_intT nRows = FS_intT(values.size()) / nPerRow;
  rows.Resize(nRows, nPerRow);
  for(FS_intT i = 0; i < nRows; ++i)
    for(FS_intT k = 0; k < nPerRow; ++k)
      rows(i, k) = renumber(values[i * nPerRow + k]);
  return rows;
}

FSIntRegisterT FSRaggedRows::Build(const FSRenumber& renumber) const
{
  FSIntRegisterT rows;
  if(counts.empty())
    return rows;

  FSIntArrayT renumbered(FS_intT(values.size()));
  for(FS_intT i = 0; i < FS_intT(values.size()); ++i)
    renumbered[i] = renumber(values[i]);
  rows.Init(ToIntArray(counts), renumbered);
  return rows;
}

void FSPolyFaceRows::Build(FSMeshPolyFaceStorage& polyFaces) const
{
  if(nFaces.empty()) {
    polyFaces.Init(0);
    polyFaces.EndNumFacesSetup();
    polyFaces.EndNumNodesSetup();
    return;
  }
  polyFaces.InitFromCountsAnd1DConnectivity(ToIntArray(nFaces), ToIntArray(nFaceNodes), ToIntArray(nodes));
}

void FSParentRows::Build(FSIntArrayT& parentOut, FSIntArrayT& parentTypeOut) const
{
  assert(parent.size() == parentType.size());
  parentOut = ToIntArray(parent);
  parentTypeOut = ToIntArray(parentType);
}

FSTopologyData FSTopologySink::Finalize() &&
{
  FSTopologyData topo;
  const FSRenumber renumber = CompactNodes(used_, nodes_, topo.globalCoords, topo.nodeGlobalNumbers);

  for(const auto& [type, rows] : inner_)
    topo.cell2NodeInner[type] = rows.Build(renumber);
  topo.cell2NodePoly2D = poly2D_.Build(renumber);
  topo.cell2NodePoly3D = poly3D_.Build(renumber);
  faces_.Build(topo.polyFaces);
  for(const auto& [type, rows] : parents_)
    rows.Build(topo.cellParent[type], topo.cellParentType[type]);
  topo.poly2DMarkers = ToIntArray(poly2DMarkers_);

  assert(topo.cell2NodePoly3D.Size() == topo.polyFaces.GetNumPolys());
  assert(topo.cell2NodePoly2D.Size() == FS_intT(poly2DMarkers_.size()));
  return topo;
}

_FS_END_NAMESPACE
