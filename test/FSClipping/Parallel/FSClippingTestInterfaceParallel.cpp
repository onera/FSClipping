#include "FSClipping/FSClippingInterfacePar.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"
#include <FSMeshCheck.h>
#include <FSMeshExportFilterHDF5.h>
#include <FSMeshExportFilterVTK.h>
#include <FSMeshExportParamsHDF5.h>
#include <FSMeshFaceExtractor.h>
#include <FSMeshImportFilterHDF5.h>
#include <FSMeshPartitionerRCB.h>
#include <FSMeshPrintInfo.h>

_FS_BEGIN_NAMESPACE

void polyMeshExportImport(FSClac* clac, FSMeshData* meshDataPtr, const FSString filename_prefix, FS_intT logLevel,
                          bool partitionIndependent)
{
  FS_intT nProcs = FSCLAC_NPROCS(clac);
  bool success;

  // export mesh as vtk
  FSMeshExportFilterVTK exportFilterVTK(clac);
  FSMeshExportParamsVTK exportParamsVTK;
  exportParamsVTK.mMeshFilename = filename_prefix + FSString(".vtk");
  exportParamsVTK.mFormatType = FSVtkEnums::FT_ASCII;
  success = exportFilterVTK.DoOp(meshDataPtr, &exportParamsVTK);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  //  --- export mesh in hdf5 format ---
  FSString filenameExport;
  FSMeshExportFilterHDF5 exportFilterHDF5(clac);
  FSMeshExportParamsHDF5 exportParamsHDF5;
  exportParamsHDF5.mExportPartitionIndependent = partitionIndependent;
  filenameExport = filename_prefix + FSString("_np");
  filenameExport.Add(nProcs);
  filenameExport.Add(".h5");
  exportParamsHDF5.mMeshFilename = filenameExport;
  success = exportFilterHDF5.DoOp(meshDataPtr, &exportParamsHDF5);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  // --- import the mesh ---
  FSMeshImportFilterHDF5 importFilterHDF5(clac);
  FSMeshImportParamsHDF5 importParamsHDF5;
  importParamsHDF5.mImportPartitionIndependent = partitionIndependent;

  importParamsHDF5.mMeshFilename = filenameExport;
  // ASSERT_TRUE(FileExists(filenameExport));
  FSMeshData meshDataImp = FSMeshData(clac);
  FSMeshData* meshDataImpPtr = &meshDataImp;
  success = importFilterHDF5.DoOp(meshDataImpPtr, &importParamsHDF5);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);
  FS_intT procID = FSCLAC_PROCID(clac);
  FSLog(clac, procID, "################## PrintInfo ###################################\n", logLevel);
  FSMeshPrintInfo printInfo(clac);
  FSMeshOpParams dummy;
  success = printInfo.DoOp(meshDataImpPtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  FSLog(clac, procID, "################## Check Mesh ###################################\n", logLevel);
  FSMeshCheck meshCheck(clac);
  success = meshCheck.DoOp(meshDataImpPtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  // --- end import mesh ---
}

static void polyMeshRepartition(FSClac* clac, FSMeshData* meshDataPtr)
{
  FSMeshOpParams dummy;
  bool success;
  FSMeshPrintInfo printInfo(clac);
  //  repartition of mesh
  FSMeshPartitionerRCB partitioner(clac);
  FSMeshPartitioningParamsRCB partitionParams;
  success = partitioner.DoOp(meshDataPtr, &partitionParams);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  // print mesh info
  success = printInfo.DoOp(meshDataPtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  // check mesh
  FSMeshCheck meshCheck(clac);
  success = meshCheck.DoOp(meshDataPtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);
}

static void polyMeshExtractFaces(FSUnstructMeshData& unstructMeshData)
{
  bool success;

  unstructMeshData.CreateLocalNumbering();

  FSMeshFaceExtractor fex;
  // don't match remote faces, keep pseudo-cells
  success = fex.PrepareFaceConnectivity(unstructMeshData, false, false);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  success = fex.PrepareFaceNodeCoordinates(FSQuantityDescArrayT());
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  FS_intT nFaces = fex.GetNLocalFaces();
  FSIntArrayT faces(nFaces, 6);
  success = fex.ExtractFaces(unstructMeshData, faces, false);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  for(FS_intT face = 0; face < nFaces; ++face) {
    const FS_intT maxFaceNodes = 7; // Maximum number of nodes that we expect of a face
    FSFloatArrayT outFaceNodeCoordinates(maxFaceNodes, 3);
    (void)fex.GetFaceNodeCoordinates(face, outFaceNodeCoordinates);
  }
}

void CheckMesh(FSClac& clac, FSMeshData* meshData)
{
  FSMeshOpParams dummy;

  FSMeshPrintInfo printInfo(&clac);
  ASSERT_TRUE(printInfo.DoOp(meshData, &dummy)) << "Print failed";

  FSMeshCheck meshCheck(&clac);
  ASSERT_TRUE(meshCheck.DoOp(meshData, &dummy)) << "Check failed";
}

TEST(FSCLippingTestInterfacePar, SurfaceInterface)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT numProc = globalClac.GetNProcs();
  FS_intT procId = globalClac.GetProcID();

  FS_intT meshID = -1;
  if(procId == 0)
    meshID = 0; // Clipper
  else if(procId % 2 == 0)
    meshID = 1;
  else
    meshID = 2;

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-8;
  FS_intT marker = -1;

  if(meshID == 1) {
    marker = 6;
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse_par.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    mesh.PrintInfo();
    polyMeshExportImport(&clac, mesh.GetMeshData(), MeshPath("output/cube_coarse_par"), 1, false);
  }

  if(meshID == 2) {
    marker = 5;
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine_par.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    mesh.PrintInfo();
    polyMeshExportImport(&clac, mesh.GetMeshData(), MeshPath("output/cube_fine_par"), 1, false);
  }

  FSClippingInterfacePar surfaceInterface(globalClac, clac, tol, marker, meshID);
  FSMesh meshClipped = surfaceInterface.BuildSurfaceInterface(mesh, meshID);

  if(meshID == 1 || meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();

    CheckMesh(clac, ptr);

    polyMeshRepartition(&clac, ptr);

    // polyMeshExtractFaces(ptr->GetUnstructCells());
  }
  if(meshID == 1) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshExportImport(&clac, ptr, MeshPath("output/cube_clipped_coarse_par"), 1, false);
  }

  if(meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshExportImport(&clac, ptr, MeshPath("output/cube_clipped_fine_par"), 1, false);
  }
}

// TEST(FSCLippingTestInterfacePar, VolumeInterface)
//{
//   FSClac globalClac(MPI_COMM_WORLD);
//   FS_intT meshID = globalClac.GetProcID();
//   FSClac clac;
//   globalClac.DivideIntoGroups(meshID, clac);
//
//   FSMesh mesh;
//   const FS_floatT tol = 1e-8;
//   FS_intT marker = -1;
//
//   if(meshID == 1) {
//     mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine.grid"));
//     marker = 1;
//   }
//
//   if(meshID == 0) {
//     mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse.grid"));
//     marker = 2;
//   }
//
//   FSClippingInterfacePar clip(globalClac, clac, tol, marker);
//   FSMesh meshClipped = clip.BuildVolumeInterface(mesh);
//
//   if(meshID == 0 || meshID == 1) {
//     FSMeshData* ptr = meshClipped.GetMeshData();
//     CheckMesh(clac, ptr);
//     polyMeshRepartition(&clac, ptr);
//   }
// }
//
//// The two cylinder boundaries only partially overlap (relative rotation):
//// faces at the rim of the interface are not fully covered by the other mesh.
//// Both sides must keep those faces whole (NOT_COVERED) so that neither
//// rebuilt mesh has holes or missing faces.
// TEST(FSCLippingTestInterfacePar, VolumeInterface2Cylinders)
//{
//   FSClac globalClac(MPI_COMM_WORLD);
//   FS_intT meshID = globalClac.GetProcID();
//   FSClac clac;
//   globalClac.DivideIntoGroups(meshID, clac);
//
//   FSMesh mesh;
//   const FS_floatT tol = 1e-6;
//   FS_intT marker = -1;
//
//   if(meshID == 0) {
//     mesh = LoadMeshWithClac(clac, MeshPath("input/mesh_cylinder_1.grid"));
//     marker = 1;
//   }
//
//   if(meshID == 1) {
//     mesh = LoadMeshWithClac(clac, MeshPath("input/mesh_cylinder_2.grid"));
//     marker = 1;
//   }
//
//   FSClippingInterfacePar clip(globalClac, clac, tol, marker);
//   FSMesh meshClipped = clip.BuildVolumeInterface(mesh);
//
//   if(meshID == 0 || meshID == 1) {
//     FSMeshData* ptr = meshClipped.GetMeshData();
//     CheckMesh(clac, ptr);
//     polyMeshRepartition(&clac, ptr);
//   }
// }

_FS_END_NAMESPACE
