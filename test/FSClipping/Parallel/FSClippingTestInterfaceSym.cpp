// Asymmetric parallel interface tests, run at exactly 6 MPI procs :
// rank 0 = clipper
// ranks 1..2 = mesh 1
// ranks 3..5 = mesh 2

#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"
#include <FSDataManagerData.h>
#include <FSDataManagerOp.h>

_FS_BEGIN_NAMESPACE

TEST(FSClippingTestInterfaceParSym, CubeVolume)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT procId = globalClac.GetProcID();
  FS_intT meshID = SymMeshID(procId);

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh(&clac);
  FSDataManagerData* data = new FSDataManagerData(&globalClac);

  FSClippedMeshParams clippedParams;
  clippedParams.mMeshKeyOriginal1 = "original1";
  clippedParams.mMeshKeyOriginal2 = "original2";
  clippedParams.mMeshKeyClipped1 = "clippedMesh1";
  clippedParams.mMeshKeyClipped2 = "clippedMesh2";
  clippedParams.mMarker1 = 6;
  clippedParams.mMarker2 = 5;
  clippedParams.mTol = 1e-8;

  if(meshID == 1) {
    FSMesh* meshPtr = data->GetMesh("original1", &clac);
    FSString meshFile = MeshPath("input/cube_hexa_coarse_par.grid");
    LoadMeshWithClac(meshPtr, meshFile);

    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, meshPtr->GetMeshData());
    mesh = *meshPtr;
  }

  if(meshID == 2) {
    FSMesh* meshPtr = data->GetMesh("original2", &clac);
    FSString meshFilename = MeshPath("input/cube_hexa_fine_par.grid");
    LoadMeshWithClac(meshPtr, meshFilename);

    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, meshPtr->GetMeshData());
    mesh = *meshPtr;
  }

  // Call the clipping operation
  FSClippedMesh clippedMesh(&globalClac);
  bool success = clippedMesh.DoOp(data, &clippedParams);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  if(meshID == 1 || meshID == 2) {
    const FSString clippedKey = meshID == 1 ? "clippedMesh1" : "clippedMesh2";
    FSMesh* meshClippedPtr = data->GetMesh(clippedKey, false);
    FSMeshData* ptr = meshClippedPtr->GetMeshData();
    if(!ptr->GetUnstructCells().HasLocalNumbering())
      ptr->GetUnstructCells().CreateLocalNumbering();
    polyMeshExtractFaces(ptr->GetUnstructCells());
    polyMeshExportImport(
      &clac, ptr, MeshPath(meshID == 1 ? "cube_clipped_vol_coarse_sym" : "output/cube_clipped_vol_fine_sym"), 1, true);
  }

  // VTK part if you want export
  // if(meshID == 1 || meshID == 2) {
  //  const FSString clippedKey = meshID == 1 ? "clippedMesh1" : "clippedMesh2";
  //  FSMesh* meshClippedPtr = data->GetMesh(clippedKey, false);
  //  if(meshClippedPtr && meshClippedPtr->IsInitialized()) {
  //    FSMeshData* ptr = meshClippedPtr->GetMeshData();
  //    // Export the clipped mesh with its as-clipped partition (one VTK per proc).
  //    ptr->GetUnstructCells().CreateLocalNumbering();
  //    exportMeshVTK(&clac, ptr,
  //                  MeshPath(meshID == 1 ? "output/cube_clipped_vol_coarse_view" :
  //                  "output/cube_clipped_vol_fine_view"));
  //  }
  //}
}

TEST(FSClippingTestInterfaceParSym, TwoCylindersVolume)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT procId = globalClac.GetProcID();
  FS_intT meshID = SymMeshID(procId);

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh(&clac);
  FSDataManagerData* data = new FSDataManagerData(&globalClac);

  FSClippedMeshParams clippedParams;
  clippedParams.mMeshKeyOriginal1 = "original1";
  clippedParams.mMeshKeyOriginal2 = "original2";
  clippedParams.mMeshKeyClipped1 = "clippedMesh1";
  clippedParams.mMeshKeyClipped2 = "clippedMesh2";
  clippedParams.mMarker1 = 1;
  clippedParams.mMarker2 = 1;
  clippedParams.mTol = 1e-6;

  if(meshID == 1) {
    FSMesh* meshPtr = data->GetMesh("original1", &clac);
    FSString meshFile = MeshPath("input/mesh_cylinder_1.grid");
    LoadMeshWithClac(meshPtr, meshFile);

    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, meshPtr->GetMeshData());
    mesh = *meshPtr;
  }

  if(meshID == 2) {
    FSMesh* meshPtr = data->GetMesh("original2", &clac);
    FSString meshFilename = MeshPath("input/mesh_cylinder_2.grid");
    LoadMeshWithClac(meshPtr, meshFilename);

    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, meshPtr->GetMeshData());
    mesh = *meshPtr;
  }

  // Call the clipping operation
  FSClippedMesh clippedMesh(&globalClac);
  bool success = clippedMesh.DoOp(data, &clippedParams);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);


  if(meshID == 1 || meshID == 2) {
    const FSString clippedKey = meshID == 1 ? "clippedMesh1" : "clippedMesh2";
    FSMesh* meshClippedPtr = data->GetMesh(clippedKey, false);
    FSMeshData* ptr = meshClippedPtr->GetMeshData();
    if(!ptr->GetUnstructCells().HasLocalNumbering())
      ptr->GetUnstructCells().CreateLocalNumbering();
    polyMeshExtractFaces(ptr->GetUnstructCells());
    polyMeshExportImport(
      &clac, ptr, MeshPath(meshID == 1 ? "output/cyl1_clipped_vol_sym" : "output/cyl2_clipped_vol_sym"), 1, true);
  }

  // VTK part if you want export
  // if(meshID == 1 || meshID == 2) {
  //  const FSString clippedKey = meshID == 1 ? "clippedMesh1" : "clippedMesh2";
  //  FSMesh* meshClippedPtr = data->GetMesh(clippedKey, false);
  //  if(meshClippedPtr && meshClippedPtr->IsInitialized()) {
  //    FSMeshData* ptr = meshClippedPtr->GetMeshData();
  //    // Export the clipped mesh with its as-clipped partition (one VTK per proc).
  //    if(!ptr->GetUnstructCells().HasLocalNumbering())
  //      ptr->GetUnstructCells().CreateLocalNumbering();
  //    exportMeshVTK(&clac, ptr,
  //                  MeshPath(meshID == 1 ? "output/cyl1_clipped_vol_view" : "output/cyl2_clipped_vol_view"));
  //  }
  //}
}
_FS_END_NAMESPACE
