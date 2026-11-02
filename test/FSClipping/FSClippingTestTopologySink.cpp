#include "gtest/gtest.h"
#include <algorithm>
#include <memory>
#include <set>
#include <vector>

#include "FSClipping/FSClippingEngine.h"
#include "FSClipping/FSTopologySink.h"
#include "FSMeshExportFilterVTK.h"
#include "FSMeshExportParamsVTK.h"
#include "MockMixedMesh.hpp"
#include "TestUtils.hpp"

_FS_BEGIN_NAMESPACE

namespace {

constexpr FS_floatT kTol = 1e-6;
constexpr FS_floatT kEps = 1e-12;

bool Near(const FSVec3& a, const FSVec3& b) { return (a - b).L2Norm() < kEps; }

// NodeId n -> newId[n]; -1 marks an unused node
FSRenumber Renumber(std::vector<FS_intT> newId) { return FSRenumber{std::move(newId)}; }

std::vector<FS_intT> RegisterRow(const FSIntRegisterT& r, FS_intT i)
{
  std::vector<FS_intT> row;
  for(FS_intT j = 0; j < r.GetNItems(i); j++)
    row.push_back(r(i, j));
  return row;
}

std::vector<FS_intT> ArrayValues(const FSIntArrayT& a)
{
  std::vector<FS_intT> v;
  for(FS_intT i = 0; i < a.Size(); i++)
    v.push_back(a[i]);
  return v;
}

void WriteVTK(FSClac& clac, FSMeshData& meshData, const std::string& filename)
{
  FSMeshExportFilterVTK exportFilter(&clac);
  FSMeshExportParamsVTK params;
  params.mMeshFilename = filename;
  params.mFormatType = FSVtkEnums::FT_ASCII;
  FSMeshData* meshDataPtr = &meshData;
  EXPECT_TRUE(exportFilter.DoOp(meshDataPtr, &params)) << "VTK export failed: " << filename;
}

// ---------------------------------------------------------------------------
// Accumulators
// ---------------------------------------------------------------------------

TEST(FSFixedRowsTest, BuildRenumbersEachRow)
{
  FSFixedRows rows{4, {0, 1, 3, 4, 4, 3, 1, 0}};
  const FSIntArrayT a = rows.Build(Renumber({0, 1, -1, 2, 3}));

  ASSERT_EQ(a.NDims(), 2);
  ASSERT_EQ(a.Size(0), 2);
  ASSERT_EQ(a.Size(1), 4);
  const std::vector<std::vector<FS_intT> > expected = {{0, 1, 2, 3}, {3, 2, 1, 0}};
  for(FS_intT i = 0; i < 2; i++)
    for(FS_intT k = 0; k < 4; k++)
      EXPECT_EQ(a(i, k), expected[i][k]) << "row " << i << ", col " << k;
}

TEST(FSFixedRowsTest, EmptyBuildsEmptyArray)
{
  FSFixedRows rows;
  EXPECT_EQ(rows.Build(Renumber({})).Size(), 0);
}

TEST(FSRaggedRowsTest, BuildKeepsCountsAndRenumbers)
{
  FSRaggedRows rows{{3, 4}, {0, 2, 5, 5, 2, 0, 6}};
  const FSIntRegisterT r = rows.Build(Renumber({0, -1, 1, -1, -1, 2, 3}));

  ASSERT_EQ(r.Size(), 2);
  EXPECT_EQ(r.ListSize(), 7);
  EXPECT_EQ(RegisterRow(r, 0), (std::vector<FS_intT>{0, 1, 2}));
  EXPECT_EQ(RegisterRow(r, 1), (std::vector<FS_intT>{2, 1, 0, 3}));
}

TEST(FSRaggedRowsTest, EmptyBuildsEmptyRegister)
{
  FSRaggedRows rows;
  EXPECT_EQ(rows.Build(Renumber({})).Size(), 0);
}

TEST(FSPolyFaceRowsTest, BuildKeepsLocalPositions)
{
  // Pyramid then tetra, faces in local positions
  FSPolyFaceRows rows;
  rows.nFaces = {5, 4};
  for(const auto& f : mock::kPyramidFaces) {
    rows.nFaceNodes.push_back(FS_intT(f.size()));
    rows.nodes.insert(rows.nodes.end(), f.begin(), f.end());
  }
  const mock::Rows tetraFaces = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
  for(const auto& f : tetraFaces) {
    rows.nFaceNodes.push_back(FS_intT(f.size()));
    rows.nodes.insert(rows.nodes.end(), f.begin(), f.end());
  }

  FSMeshPolyFaceStorage storage;
  rows.Build(storage);

  ASSERT_EQ(storage.GetNumPolys(), 2);
  const std::vector<const mock::Rows*> expected = {&mock::kPyramidFaces, &tetraFaces};
  for(FS_intT p = 0; p < 2; p++) {
    const mock::Rows& faces = *expected[p];
    ASSERT_EQ(storage.GetNFaces(p), FS_intT(faces.size()));
    for(FS_intT f = 0; f < FS_intT(faces.size()); f++) {
      ASSERT_EQ(storage.GetNFaceNodes(p, f), FS_intT(faces[f].size()));
      for(FS_intT k = 0; k < FS_intT(faces[f].size()); k++)
        EXPECT_EQ(storage.GetFaceNode(p, f, k), faces[f][k]) << "poly " << p << ", face " << f;
    }
  }
}

TEST(FSPolyFaceRowsTest, EmptyLeavesStorageEmpty)
{
  FSPolyFaceRows rows;
  FSMeshPolyFaceStorage storage;
  rows.Build(storage);
  EXPECT_EQ(storage.GetNumPolys(), 0);
}

TEST(FSParentRowsTest, BuildCopiesParentsAndTypes)
{
  FSParentRows rows{{12, 7, 3}, {FSMeshEnums::CT_Hexa8, FSMeshEnums::CT_Tetra4, FSMeshEnums::CT_Poly3D}};
  FSIntArrayT parent, parentType;
  rows.Build(parent, parentType);

  EXPECT_EQ(ArrayValues(parent), rows.parent);
  EXPECT_EQ(ArrayValues(parentType), rows.parentType);
}

// ---------------------------------------------------------------------------
// Sink on the mixed mock mesh: hexa, Poly3D pyramid, tetra, Quad4, Tri3, Poly2D
// ---------------------------------------------------------------------------

class FSTopologySinkTest : public ::testing::Test
{
protected:
  FSClac clac;
  FSMesh mesh = mock::BuildMixedMesh(clac);
  std::vector<FSCellSource> sources = BuildCellSources(mesh);
  std::vector<FSFaceMatch> matches;
  BoundaryExtraction be;
  std::unique_ptr<FSCellClipped> clipped;
  std::unique_ptr<FSNodeRegistry> nodes;

  const FSCellSource& Source(FSMeshEnums::CellType t) const
  {
    return *std::find_if(sources.begin(), sources.end(), [t](const FSCellSource& s) { return s.Type() == t; });
  }

  FS_intT Cell(FSMeshEnums::CellType t) const { return mock::Offset(mesh, t); }

  FS_intT LocalFace(FSMeshEnums::CellType t, const std::set<FS_intT>& corners) const
  {
    const FSCellSource& src = Source(t);
    const FS_intT c = Cell(t);
    for(FS_intT f = 0; f < src.NFaces(c); f++) {
      std::set<FS_intT> s;
      for(FS_intT k = 0; k < src.NFaceCorners(c, f); k++)
        s.insert(src.FaceCorner(c, f, k));
      if(s == corners)
        return f;
    }
    ADD_FAILURE() << "no such face";
    return -1;
  }

  // Hexa bottom face {0,1,2,3} split in two by P(.5,0,0) - R(.5,1,0): P and R are new NodeIds 10 and 11
  void SplitHexaBottom()
  {
    const auto t = FSMeshEnums::CT_Hexa8;
    matches.push_back(
      mock::MakeMatch(Cell(t), t, {{0, 0, 0}, {0.5, 0, 0}, {0.5, 1, 0}, {0, 1, 0}}, {10, 100, 101, 13}));
    matches.push_back(
      mock::MakeMatch(Cell(t), t, {{0.5, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0.5, 1, 0}}, {100, 11, 12, 101}));
    be.cell2BoundaryFaces[Cell(t)].push_back(LocalFace(t, {0, 1, 2, 3}));
    // The Quad4 closing that face is replaced by the match polygons, as FSBoundaryFaceProvider records it
    be.surfaceCells.insert(Cell(FSMeshEnums::CT_Quad4));
  }

  // Fields the clipper fills on each match: closed surface cell and its carried marker
  void FillClipperFields(FS_intT marker)
  {
    for(size_t m = 0; m < matches.size(); m++) {
      matches[m].faceOwner1 = Cell(FSMeshEnums::CT_Quad4);
      matches[m].faceOwnerType1 = FSMeshEnums::CT_Quad4;
      matches[m].faceMarker1 = marker;
    }
  }

  void BuildRegistry()
  {
    clipped = std::make_unique<FSCellClipped>(matches, be);
    nodes = std::make_unique<FSNodeRegistry>(mesh, sources, *clipped, kTol);
  }

  // Coordinates behind each entry of a row of the output
  std::vector<FSVec3> RowCoords(const FSTopologyData& topo, const std::vector<FS_intT>& row) const
  {
    std::vector<FSVec3> c;
    for(FS_intT n : row)
      c.push_back(topo.globalCoords[n]);
    return c;
  }

  std::vector<FSVec3> RegistryCoords(const std::vector<NodeId>& ids) const
  {
    std::vector<FSVec3> c;
    for(NodeId n : ids)
      c.push_back(nodes->Coord(n));
    return c;
  }

  static void ExpectNear(const std::vector<FSVec3>& a, const std::vector<FSVec3>& b)
  {
    ASSERT_EQ(a.size(), b.size());
    for(size_t i = 0; i < a.size(); i++)
      EXPECT_TRUE(Near(a[i], b[i])) << "entry " << i;
  }
};

TEST_F(FSTopologySinkTest, CompactNodesKeepsUsedNodesInOrder)
{
  SplitHexaBottom();
  BuildRegistry();
  ASSERT_EQ(nodes->Size(), 12);

  std::vector<bool> used(12, false);
  for(NodeId n : {1, 3, 4, 10})
    used[n] = true;

  std::vector<FSVec3> coords;
  FSIntArrayT clipperIds;
  const FSRenumber r = CompactNodes(used, *nodes, coords, clipperIds);

  EXPECT_EQ(r.newId, (std::vector<FS_intT>{-1, 0, -1, 1, 2, -1, -1, -1, -1, -1, 3, -1}));
  ExpectNear(coords, RegistryCoords({1, 3, 4, 10}));
  EXPECT_EQ(ArrayValues(clipperIds), (std::vector<FS_intT>{11, 13, -1, 100}));
}

TEST_F(FSTopologySinkTest, CompactNodesNothingUsed)
{
  BuildRegistry();
  std::vector<FSVec3> coords;
  FSIntArrayT clipperIds;
  const FSRenumber r = CompactNodes(std::vector<bool>(nodes->Size(), false), *nodes, coords, clipperIds);

  EXPECT_TRUE(coords.empty());
  EXPECT_EQ(clipperIds.Size(), 0);
  EXPECT_TRUE(std::all_of(r.newId.begin(), r.newId.end(), [](FS_intT k) { return k == -1; }));
}

TEST_F(FSTopologySinkTest, AddStandardDropsUnusedNodes)
{
  BuildRegistry();
  FSTopologySink sink(*nodes);
  sink.AddStandard(Source(FSMeshEnums::CT_Tetra4), Cell(FSMeshEnums::CT_Tetra4));
  FSTopologyData topo = std::move(sink).Finalize();

  // Only the tetra nodes {5, 6, 8, 9} are kept, in registry order
  ExpectNear(topo.globalCoords, RegistryCoords({5, 6, 8, 9}));
  EXPECT_EQ(ArrayValues(topo.nodeGlobalNumbers), (std::vector<FS_intT>{-1, -1, -1, -1}));

  const FSIntArrayT& tetra = topo.cell2NodeInner[FSMeshEnums::CT_Tetra4];
  ASSERT_EQ(tetra.Size(0), 1);
  ASSERT_EQ(tetra.Size(1), 4);
  for(FS_intT k = 0; k < 4; k++)
    EXPECT_EQ(tetra(0, k), k);

  EXPECT_EQ(ArrayValues(topo.cellParent[FSMeshEnums::CT_Tetra4]), (std::vector<FS_intT>{Cell(FSMeshEnums::CT_Tetra4)}));
  EXPECT_EQ(ArrayValues(topo.cellParentType[FSMeshEnums::CT_Tetra4]), (std::vector<FS_intT>{FSMeshEnums::CT_Tetra4}));
  EXPECT_EQ(topo.cell2NodePoly2D.Size(), 0);
  EXPECT_EQ(topo.cell2NodePoly3D.Size(), 0);
  EXPECT_EQ(topo.polyFaces.GetNumPolys(), 0);
}

TEST_F(FSTopologySinkTest, AddPoly3DFromSourceCopiesRowAndFaces)
{
  BuildRegistry();
  FSTopologySink sink(*nodes);
  const FS_intT pyramid = Cell(FSMeshEnums::CT_Poly3D);
  sink.AddPoly3D(Source(FSMeshEnums::CT_Poly3D), pyramid);
  FSTopologyData topo = std::move(sink).Finalize();

  ASSERT_EQ(topo.cell2NodePoly3D.Size(), 1);
  ExpectNear(RowCoords(topo, RegisterRow(topo.cell2NodePoly3D, 0)), RegistryCoords(mock::kPyramid));

  // Local positions are the original ones
  ASSERT_EQ(topo.polyFaces.GetNumPolys(), 1);
  ASSERT_EQ(topo.polyFaces.GetNFaces(0), FS_intT(mock::kPyramidFaces.size()));
  for(FS_intT f = 0; f < FS_intT(mock::kPyramidFaces.size()); f++)
    for(FS_intT k = 0; k < FS_intT(mock::kPyramidFaces[f].size()); k++)
      EXPECT_EQ(topo.polyFaces.GetFaceNode(0, f, k), mock::kPyramidFaces[f][k]);

  EXPECT_EQ(ArrayValues(topo.cellParent[FSMeshEnums::CT_Poly3D]), (std::vector<FS_intT>{pyramid}));
  EXPECT_EQ(ArrayValues(topo.cellParentType[FSMeshEnums::CT_Poly3D]), (std::vector<FS_intT>{FSMeshEnums::CT_Poly3D}));
}

TEST_F(FSTopologySinkTest, AddPoly3DFromPolyCellMapsFacesThroughTheRow)
{
  SplitHexaBottom();
  BuildRegistry();
  const auto t = FSMeshEnums::CT_Hexa8;
  const FSPolyCell poly = BuildClippedCell(Source(t), Cell(t), *clipped, *nodes);
  ASSERT_EQ(poly.faces.size(), 7u);

  FSTopologySink sink(*nodes);
  sink.AddPoly3D(Cell(t), t, poly);
  FSTopologyData topo = std::move(sink).Finalize();

  // 8 hexa nodes + P + R, each node once in the row
  ASSERT_EQ(topo.cell2NodePoly3D.Size(), 1);
  const std::vector<FS_intT> row = RegisterRow(topo.cell2NodePoly3D, 0);
  EXPECT_EQ(row.size(), 10u);
  EXPECT_EQ(std::set<FS_intT>(row.begin(), row.end()).size(), 10u);
  EXPECT_EQ(topo.globalCoords.size(), 10u);

  // Each face, read through polyFaces and the row, gives back the FSPolyCell face
  ASSERT_EQ(topo.polyFaces.GetNumPolys(), 1);
  ASSERT_EQ(topo.polyFaces.GetNFaces(0), 7);
  for(FS_intT f = 0; f < 7; f++) {
    std::vector<FS_intT> face;
    for(FS_intT k = 0; k < topo.polyFaces.GetNFaceNodes(0, f); k++)
      face.push_back(row[topo.polyFaces.GetFaceNode(0, f, k)]);
    ExpectNear(RowCoords(topo, face), RegistryCoords(poly.faces[f]));
  }

  EXPECT_EQ(ArrayValues(topo.cellParent[FSMeshEnums::CT_Poly3D]), (std::vector<FS_intT>{Cell(t)}));
  EXPECT_EQ(ArrayValues(topo.cellParentType[FSMeshEnums::CT_Poly3D]), (std::vector<FS_intT>{t}));
}

TEST_F(FSTopologySinkTest, AddPoly2DStoresRowsParentsAndMarkers)
{
  SplitHexaBottom();
  BuildRegistry();
  const FS_intT quad = Cell(FSMeshEnums::CT_Quad4);

  FSTopologySink sink(*nodes);
  clipped->ForEachMatch([&](FS_intT m, const FSFaceMatch&) {
    sink.AddPoly2D(quad, FSMeshEnums::CT_Quad4, nodes->Polygon(m), 70 + m);
  });
  FSTopologyData topo = std::move(sink).Finalize();

  // Nodes 0..3 of the bottom face, then P and R
  ExpectNear(topo.globalCoords, RegistryCoords({0, 1, 2, 3, 10, 11}));
  EXPECT_EQ(ArrayValues(topo.nodeGlobalNumbers), (std::vector<FS_intT>{10, 11, 12, 13, 100, 101}));

  ASSERT_EQ(topo.cell2NodePoly2D.Size(), 2);
  for(FS_intT m = 0; m < 2; m++) {
    const auto polygon = nodes->Polygon(m);
    ExpectNear(RowCoords(topo, RegisterRow(topo.cell2NodePoly2D, m)),
               RegistryCoords(std::vector<NodeId>(polygon.begin(), polygon.end())));
  }
  EXPECT_EQ(ArrayValues(topo.poly2DMarkers), (std::vector<FS_intT>{70, 71}));
  EXPECT_EQ(ArrayValues(topo.cellParent[FSMeshEnums::CT_Poly2D]), (std::vector<FS_intT>{quad, quad}));
  EXPECT_EQ(ArrayValues(topo.cellParentType[FSMeshEnums::CT_Poly2D]),
            (std::vector<FS_intT>{FSMeshEnums::CT_Quad4, FSMeshEnums::CT_Quad4}));
}

// Every cell of the mixed mesh, emitted as the single loop of BuildTopology will
TEST_F(FSTopologySinkTest, MixedMeshKeepsContainersAligned)
{
  SplitHexaBottom();
  BuildRegistry();
  FSTopologySink sink(*nodes);

  for(const auto& src : sources)
    src.ForEachOwnedCell([&](FS_intT c) {
      if(clipped->HasBoundaryFace(c))
        sink.AddPoly3D(c, src.Type(), BuildClippedCell(src, c, *clipped, *nodes));
      else if(src.Type() == FSMeshEnums::CT_Quad4)
        return; // hexa bottom: replaced by its polygons
      else if(src.IsPoly3D())
        sink.AddPoly3D(src, c);
      else if(src.IsPoly2D())
        sink.AddPoly2D(c, src.Type(), src.Row(c), -1);
      else
        sink.AddStandard(src, c);
    });
  clipped->ForEachMatch([&](FS_intT m, const FSFaceMatch&) {
    sink.AddPoly2D(Cell(FSMeshEnums::CT_Quad4), FSMeshEnums::CT_Quad4, nodes->Polygon(m), 70 + m);
  });
  FSTopologyData topo = std::move(sink).Finalize();

  // All 10 old nodes plus P and R are used
  EXPECT_EQ(topo.globalCoords.size(), 12u);
  EXPECT_EQ(topo.nodeGlobalNumbers.Size(), 12);

  // Clipped hexa and copied pyramid
  EXPECT_EQ(topo.cell2NodePoly3D.Size(), 2);
  EXPECT_EQ(topo.polyFaces.GetNumPolys(), 2);
  EXPECT_EQ(ArrayValues(topo.cellParentType[FSMeshEnums::CT_Poly3D]),
            (std::vector<FS_intT>{FSMeshEnums::CT_Hexa8, FSMeshEnums::CT_Poly3D}));

  // Original Poly2D first, then the two match polygons
  EXPECT_EQ(topo.cell2NodePoly2D.Size(), 3);
  EXPECT_EQ(ArrayValues(topo.poly2DMarkers), (std::vector<FS_intT>{-1, 70, 71}));
  EXPECT_EQ(ArrayValues(topo.cellParentType[FSMeshEnums::CT_Poly2D]),
            (std::vector<FS_intT>{FSMeshEnums::CT_Poly2D, FSMeshEnums::CT_Quad4, FSMeshEnums::CT_Quad4}));

  // Tetra and Tri3 stay standard; the replaced Quad4 is gone
  EXPECT_EQ(topo.cell2NodeInner[FSMeshEnums::CT_Tetra4].Size(0), 1);
  EXPECT_EQ(topo.cell2NodeInner[FSMeshEnums::CT_Tri3].Size(0), 1);
  EXPECT_FALSE(topo.cell2NodeInner.Has(FSMeshEnums::CT_Quad4));
  EXPECT_FALSE(topo.cell2NodeInner.Has(FSMeshEnums::CT_Hexa8));
}

// Every cell of the mixed mesh, emitted as the single loop of BuildTopology will
TEST_F(FSTopologySinkTest, FSClippingEngine)
{

  FS_floatT tol = 1e-8;
  SplitHexaBottom();
  FSClippingEngine engine(tol);

  FSTopologyData topo = engine.BuildTopology(mesh, be, matches);
  // All 10 old nodes plus P and R are used

  auto meshClipped = engine.Reconstruct(clac, mesh, topo);
  meshClipped.Check();

  WriteVTK(clac, *mesh.GetMeshData(), MeshPath("output/engine_mock_orig.vtk"));
  WriteVTK(clac, *meshClipped.GetMeshData(), MeshPath("output/engine_mock_clipped.vtk"));
}

// CADGroupID: kept cells read their parent, match polygons take the marker carried by the clipper
TEST_F(FSTopologySinkTest, Poly2DMarkersAreCopied)
{
  const FS_floatT tol = 1e-8;
  const FSString cadGroupID = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  auto addMarker = [&](FSMeshEnums::CellType t, FS_intT v) {
    FSIntArrayT values(1);
    values[0] = v;
    mesh.GetMeshData()->GetUnstructCells().InitCellAttribute(cadGroupID, t, values);
  };
  addMarker(FSMeshEnums::CT_Quad4, 1);
  addMarker(FSMeshEnums::CT_Tri3, 2);
  addMarker(FSMeshEnums::CT_Poly2D, 3);

  // Carried marker differs from the Quad4 one, to tell both paths apart
  SplitHexaBottom();
  FillClipperFields(7);

  FSClippingEngine engine(tol);
  FSTopologyData topo = engine.BuildTopology(mesh, be, matches);
  EXPECT_EQ(ArrayValues(topo.poly2DMarkers), (std::vector<FS_intT>{-1, 7, 7}));
  FSMesh meshClipped = engine.Reconstruct(clac, mesh, topo);

  const FSUnstructMeshData& um = meshClipped.GetMeshData()->GetUnstructCells();
  ASSERT_TRUE(um.HasCellAttribute(cadGroupID, FSMeshEnums::CT_Poly2D)) << "CADGroupID not copied to Poly2D";
  EXPECT_EQ(ArrayValues(um.GetCellAttribute(cadGroupID, FSMeshEnums::CT_Poly2D)), (std::vector<FS_intT>{3, 7, 7}));
  ASSERT_TRUE(um.HasCellAttribute(cadGroupID, FSMeshEnums::CT_Tri3)) << "CADGroupID not copied to Tri3";
  EXPECT_EQ(ArrayValues(um.GetCellAttribute(cadGroupID, FSMeshEnums::CT_Tri3)), (std::vector<FS_intT>{2}));
}

} // namespace

_FS_END_NAMESPACE
