/*
Attributes of the clipped Poly2D cells, run at exactly 9 MPI procs. It is the smallest count at which RCB splits a
volume/surface pair on these meshes, which is what this test is about.

A Poly2D is born from a match, and its parent is the surface cell the matched face closes. When RCB splits a
volume/surface pair, that parent is owned by another proc: its index means nothing in the local pool, so the
CADGroupID cannot be read back locally and has to travel with the match. Without that, the copy either reads an
unrelated local cell or runs past the end of the pool.

At 9 procs, the multi-element mesh also leaves a proc without any Poly2D: its empty pool must still carry the same
attributes as the others, or the parallel HDF5 export aborts.
*/

#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"
#include <FSDataManagerData.h>
#include <FSMeshImportFilterHDF5.h>

_FS_BEGIN_NAMESPACE

namespace {

/// Counts the Poly2D cells of `mesh` whose CADGroupID differs from `expectedMarker`.
FS_int32T CountWrongPoly2DMarkers(const FSMesh& mesh, FS_intT expectedMarker)
{
  static const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  const FSUnstructMeshData& meshData = mesh.GetMeshData()->GetUnstructCells();

  const FSCellPool* cellPool = meshData.GetCellPool(FSMeshEnums::CT_Poly2D);
  if(!cellPool || meshData.GetCellAttributes(FSMeshEnums::CT_Poly2D).Size() == 0)
    return 0;

  const FSIntArrayT& markers = meshData.GetCellAttribute(attrCADGroupIDName, FSMeshEnums::CT_Poly2D);

  FS_int32T wrong = 0;
  for(FS_intT i = 0; i < markers.Size(); ++i)
    if(markers[i] != expectedMarker)
      ++wrong;
  return wrong;
}

/// Returns 1 when the Poly2D pool of `mesh` lacks the CADGroupID attribute, 0 otherwise.
FS_int32T CountMissingPoly2DMarker(const FSMesh& mesh)
{
  static const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  return mesh.GetMeshData()->GetUnstructCells().HasCellAttribute(attrCADGroupIDName, FSMeshEnums::CT_Poly2D) ? 0 : 1;
}

/// One input mesh of the clipping: file and boundary marker.
struct FSClippingInput {
  const char* path;
  FS_intT marker;
};

class FSClippingTestPoly2DAttributesPar : public ::testing::Test
{
protected:
  /// Clips `input1` against `input2`; returns false when the proc count is not 9.
  bool RunClipping(FSClippingInput input1, FSClippingInput input2)
  {
    if(globalClac.GetNProcs() != 9)
      return false;

    meshID = SymMeshID(globalClac.GetProcID());
    globalClac.DivideIntoGroups(meshID, clac);
    data = new FSDataManagerData(&globalClac);

    params.mMeshKeyOriginal1 = "original1";
    params.mMeshKeyOriginal2 = "original2";
    params.mMeshKeyClipped1 = "clippedMesh1";
    params.mMeshKeyClipped2 = "clippedMesh2";
    params.mMarker1 = input1.marker;
    params.mMarker2 = input2.marker;
    params.mTol = 1e-8;

    if(meshID == 1 || meshID == 2) {
      FSMesh* meshPtr = data->GetMesh(meshID == 1 ? "original1" : "original2", &clac);
      const FSString path = MeshPath((meshID == 1 ? input1 : input2).path);
      if(path.EndsWith(".h5")) {
        FSMeshImportParamsHDF5 importParams;
        importParams.mMeshFilename = path;
        EXPECT_TRUE(meshPtr->ImportMesh(&importParams));
      } else
        LoadMeshWithClac(meshPtr, path);
      polyMeshRepartition(&clac, meshPtr->GetMeshData());
      meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    }

    FSClippedMesh clippedMesh(&globalClac);
    EXPECT_TRUE(clippedMesh.DoOp(data, &params));
    return true;
  }

  const FSMesh& Clipped() { return *data->GetMesh(meshID == 1 ? "clippedMesh1" : "clippedMesh2", &clac); }

  FSClac globalClac{MPI_COMM_WORLD};
  FSClac clac;
  FS_intT meshID = -1;
  FSDataManagerData* data = nullptr;
  FSClippedMeshParams params;
};

constexpr FSClippingInput kCoarseCube{"input/cube_hexa_coarse_par.grid", 6};
constexpr FSClippingInput kFineCube{"input/cube_hexa_fine_par.grid", 5};
// Hexa + Poly3D + Poly2D: at 9 procs, some proc of this mesh owns no Poly2D
constexpr FSClippingInput kMultiElement{"input/cube_hexa_translate.h5", 1};
constexpr FSClippingInput kFineCubeSerial{"input/cube_hexa_fine.grid", 2};

} // namespace

TEST_F(FSClippingTestPoly2DAttributesPar, ClippedFacesKeepTheirCADGroupID)
{
  if(!RunClipping(kCoarseCube, kFineCube))
    GTEST_SKIP() << "requires exactly 9 MPI procs";
  if(meshID == 0)
    return; // clipper proc holds no clipped mesh

  const FS_intT expectedMarker = (meshID == 1) ? params.mMarker1 : params.mMarker2;
  FS_int32T localWrong = CountWrongPoly2DMarkers(Clipped(), expectedMarker);
  FS_int32T totalWrong = 0;
  clac.Sum(&localWrong, &totalWrong, 1);

  // Every clipped face lies on the clipped boundary, so all Poly2D must carry that marker
  EXPECT_EQ(totalWrong, 0);
}

TEST_F(FSClippingTestPoly2DAttributesPar, EveryProcCarriesThePoly2DMarker)
{
  if(!RunClipping(kMultiElement, kFineCubeSerial))
    GTEST_SKIP() << "requires exactly 9 MPI procs";
  if(meshID == 0)
    return;

  // A proc with no Poly2D still has the CADGroupID attribute, empty, on its Poly2D pool
  FS_int32T localMissing = CountMissingPoly2DMarker(Clipped());
  FS_int32T totalMissing = 0;
  clac.Sum(&localMissing, &totalMissing, 1);
  EXPECT_EQ(totalMissing, 0);
}

TEST_F(FSClippingTestPoly2DAttributesPar, ExportImportHDF5WithAnEmptyPoly2DPool)
{
  if(!RunClipping(kMultiElement, kFineCubeSerial))
    GTEST_SKIP() << "requires exactly 9 MPI procs";
  if(meshID == 0)
    return;

  // The CADGroupID missing on the empty pool made the parallel HDF5 export abort in MPI_Unpack
  FSMesh& clipped = *data->GetMesh(meshID == 1 ? "clippedMesh1" : "clippedMesh2", &clac);
  polyMeshExportImport(&clac, clipped.GetMeshData(),
                       MeshPath(meshID == 1 ? "output/poly2DAttributes1" : "output/poly2DAttributes2"), 1, false);
}

_FS_END_NAMESPACE
