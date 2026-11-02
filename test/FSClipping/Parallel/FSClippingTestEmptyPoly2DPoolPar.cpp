/*
Clipped mesh whose procs do not hold the same cell types, run at exactly 4 MPI procs: 1 clipper, cube_hexa_coarse on
2 procs, cube_hexa_fine on 1 proc. With the default RCB, one proc of the coarse cube owns no Poly2D and no Hexa8, only
boundary Quad4. The reconstruction must still give it the same pools and attributes as the other proc: a missing
Hexa8 pool shifts the collective pool inits (deadlock), a missing CADGroupID makes the parallel HDF5 export abort.
*/

#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"
#include <FSDataManagerData.h>

_FS_BEGIN_NAMESPACE

namespace {

const FSMeshEnums::CellType kCheckedTypes[] = {FSMeshEnums::CT_Quad4, FSMeshEnums::CT_Hexa8, FSMeshEnums::CT_Poly2D,
                                               FSMeshEnums::CT_Poly3D};

FS_int32T NCells(const FSUnstructMeshData& um, FSMeshEnums::CellType t)
{
  return um.GetCellPool(t) ? static_cast<FS_int32T>(um.GetNCells(t)) : 0;
}

class FSClippingTestEmptyPoly2DPoolPar : public ::testing::Test
{
protected:
  /// Clips the coarse cube against the fine one; returns false when the proc count is not 4.
  bool RunClipping()
  {
    if(globalClac.GetNProcs() != 4)
      return false;

    meshID = AsymMeshID(globalClac.GetProcID(), 2);
    globalClac.DivideIntoGroups(meshID, clac);
    data = new FSDataManagerData(&globalClac);

    FSClippedMeshParams params;
    params.mMeshKeyOriginal1 = "original1";
    params.mMeshKeyOriginal2 = "original2";
    params.mMeshKeyClipped1 = "clippedMesh1";
    params.mMeshKeyClipped2 = "clippedMesh2";
    params.mMarker1 = 2;
    params.mMarker2 = 1;
    params.mTol = 1e-8;

    if(meshID == 1 || meshID == 2) {
      FSMesh* meshPtr = data->GetMesh(meshID == 1 ? "original1" : "original2", &clac);
      LoadMeshWithClac(meshPtr, MeshPath(meshID == 1 ? "input/cube_hexa_coarse.grid" : "input/cube_hexa_fine.grid"));
      polyMeshRepartition(&clac, meshPtr->GetMeshData());
      meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    }

    FSClippedMesh clippedMesh(&globalClac);
    EXPECT_TRUE(clippedMesh.DoOp(data, &params));
    return true;
  }

  FSUnstructMeshData& Clipped1() { return data->GetMesh("clippedMesh1", &clac)->GetMeshData()->GetUnstructCells(); }

  /// Fails unless one coarse proc owns no Poly2D and no Hexa8 while the other owns some: otherwise nothing is tested.
  void ExpectAsymmetricSplit()
  {
    const FSUnstructMeshData& um = Clipped1();
    for(const FSMeshEnums::CellType t : {FSMeshEnums::CT_Poly2D, FSMeshEnums::CT_Hexa8}) {
      FS_int32T n = NCells(um, t), nMin = 0, nMax = 0;
      clac.GetMin(&n, &nMin, 1);
      clac.GetMax(&n, &nMax, 1);
      EXPECT_EQ(nMin, 0) << FSMeshEnums::CellTypeToString(t) << " on every coarse proc";
      EXPECT_GT(nMax, 0) << FSMeshEnums::CellTypeToString(t) << " on no coarse proc";
    }
  }

  FSClac globalClac{MPI_COMM_WORLD};
  FSClac clac;
  FS_intT meshID = -1;
  FSDataManagerData* data = nullptr;
};

} // namespace

TEST_F(FSClippingTestEmptyPoly2DPoolPar, EveryProcHasTheSamePoolsAndAttributes)
{
  if(!RunClipping())
    GTEST_SKIP() << "requires exactly 4 MPI procs";
  if(meshID != 1)
    return;
  ExpectAsymmetricSplit();

  const FSUnstructMeshData& um = Clipped1();
  for(const FSMeshEnums::CellType t : kCheckedTypes) {
    FS_int32T hasPool = um.GetCellPool(t) ? 1 : 0, nWithPool = 0;
    clac.Sum(&hasPool, &nWithPool, 1);
    EXPECT_EQ(nWithPool, clac.GetNProcs()) << FSMeshEnums::CellTypeToString(t) << " pool missing on a proc";

    FS_int32T nAttribs = hasPool ? static_cast<FS_int32T>(um.GetCellAttributes(t).Size()) : 0, nMin = 0, nMax = 0;
    clac.GetMin(&nAttribs, &nMin, 1);
    clac.GetMax(&nAttribs, &nMax, 1);
    EXPECT_EQ(nMin, nMax) << FSMeshEnums::CellTypeToString(t) << " attributes differ between procs";
  }

  static const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  FS_int32T hasMarker = um.HasCellAttribute(attrCADGroupIDName, FSMeshEnums::CT_Poly2D) ? 1 : 0, nWithMarker = 0;
  clac.Sum(&hasMarker, &nWithMarker, 1);
  EXPECT_EQ(nWithMarker, clac.GetNProcs());
}

TEST_F(FSClippingTestEmptyPoly2DPoolPar, ExportImportHDF5)
{
  if(!RunClipping())
    GTEST_SKIP() << "requires exactly 4 MPI procs";
  if(meshID != 1)
    return;
  ExpectAsymmetricSplit();

  // An empty pool lacking CADGroupID made the parallel HDF5 export abort in MPI_Unpack
  polyMeshExportImport(&clac, data->GetMesh("clippedMesh1", &clac)->GetMeshData(), MeshPath("output/emptyPoly2DPool"),
                       1, false);
}

_FS_END_NAMESPACE
