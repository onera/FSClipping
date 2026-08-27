#include "FSMeshImportParamsTAU.h"
#include "FSMeshExportFilterVTK.h"
#include "FSMeshExportFilterVTK.h"
#include "FSMeshExportParamsVTK.h"
#include "FSMeshExportFilterHDF5.h"
#include "FSMeshImportFilterHDF5.h"
#include "FSMeshPartitionerRCB.h"
#include "gtest/gtest.h"
#include <FSMeshCheck.h>
#include <FSMeshData.h>
#include <FSMeshPrintInfo.h>

#include "FSClipping/FSClippingInterface.h"
#include "TestUtils.hpp"

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

namespace {
FSMesh LoadMeshTau(FSClac& clac, const std::string& filename)
{
  FSMeshImportParamsTAU params;
  params.mMeshFilename = filename;

  FSMesh mesh(&clac);
  EXPECT_TRUE(mesh.ImportMesh(&params));

  mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  return mesh;
}

FSMesh LoadMeshHDF5(FSClac& clac, const std::string& filename)
{
  FSMeshImportParamsHDF5 params;
  params.mMeshFilename = filename;

  FSMesh mesh(&clac);
  EXPECT_TRUE(mesh.ImportMesh(&params));

  mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  return mesh;
}

void CheckMesh(FSClac& clac, FSMeshData* meshData)
{
  FSMeshOpParams dummy;

  FSMeshPrintInfo printInfo(&clac);
  ASSERT_TRUE(printInfo.DoOp(meshData, &dummy)) << "Print failed";

  FSMeshCheck meshCheck(&clac);
  ASSERT_TRUE(meshCheck.DoOp(meshData, &dummy)) << "Check failed";
}
} // namespace

TEST(FSClippingTestInterface, ReconstructClippedMesh2DHexaHexa)
{
  FSClac clac1, clac2, clac3;

  auto mesh1 = LoadMeshTau(clac1, MeshPath("input/cube_hexa_coarse.grid"));
  auto mesh2 = LoadMeshTau(clac2, MeshPath("input/cube_hexa_fine.grid"));

  const FS_intT marker1 = 2;
  const FS_intT marker2 = 1;
  const FS_floatT tol = 1e-8;

  FSClippingInterface clip(clac1, clac2, tol, marker1, marker2);

  auto meshClipped = clip.BuildSurfaceInterface(mesh1, mesh2);
  FSMeshData* ptr = meshClipped.GetMeshData();

  CheckMesh(clac3, ptr);

  polyMeshRepartition(&clac3, ptr);

  polyMeshExtractFaces(ptr->GetUnstructCells());

  polyMeshExportImport(&clac3, ptr, MeshPath("output/cube2D_clipped_HexaHexa"), 1, false);
}

TEST(FSClippingTestInterface, ReconstructClippedMesh3DHexaHexa)
{
  FSClac clac1, clac2, clac3;

  auto mesh1 = LoadMeshTau(clac1, MeshPath("input/cube_hexa_coarse.grid"));
  auto mesh2 = LoadMeshTau(clac2, MeshPath("input/cube_hexa_fine.grid"));

  const FS_intT marker1 = 2;
  const FS_intT marker2 = 1;
  const FS_floatT tol = 1e-8;

  FSClippingInterface clip(clac1, clac2, tol, marker1, marker2);

  auto meshClipped = clip.BuildVolumeInterface(mesh1, mesh2);
  FSMeshData* ptr = meshClipped.GetMeshData();

  CheckMesh(clac3, ptr);

  polyMeshRepartition(&clac3, ptr);

  polyMeshExtractFaces(ptr->GetUnstructCells());

  polyMeshExportImport(&clac3, ptr, MeshPath("output/cube3D_clipped_HexaHexa"), 1, false);
}

TEST(FSClippingTestInterface, ReconstructClippedMesh3DHexaTetra)
{
  FSClac clac1, clac2, clac3;

  auto mesh1 = LoadMeshTau(clac1, MeshPath("input/cube_hexa_coarse.grid"));
  auto mesh2 = LoadMeshTau(clac2, MeshPath("input/cube_tetra_fine.grid"));

  const FS_intT marker1 = 2;
  const FS_intT marker2 = 1;
  const FS_floatT tol = 1e-8; // hole when tol < 1e-3

  FSClippingInterface clip(clac1, clac2, tol, marker1, marker2);

  auto meshClipped = clip.BuildVolumeInterface(mesh1, mesh2);
  FSMeshData* ptr = meshClipped.GetMeshData();

  CheckMesh(clac3, ptr);

  polyMeshRepartition(&clac3, ptr);

  polyMeshExtractFaces(ptr->GetUnstructCells());

  polyMeshExportImport(&clac3, ptr, MeshPath("output/cube3D_clipped_HexaTetra"), 1, false);
}

TEST(FSClippingTestInterface, ReconstructClippedMesh3DHexaHexaRotate)
{
  FSClac clac1, clac2, clac3;

  auto mesh1 = LoadMeshTau(clac1, MeshPath("input/square_mesh.grid"));
  auto mesh2 = LoadMeshTau(clac2, MeshPath("input/square_mesh_rotate.grid"));

  const FS_intT marker1 = 2;
  const FS_intT marker2 = 1;
  const FS_floatT tol = 1e-10;

  FSClippingInterface clip(clac1, clac2, tol, marker1, marker2);

  auto meshClipped = clip.BuildVolumeInterface(mesh1, mesh2);
  FSMeshData* ptr = meshClipped.GetMeshData();

  CheckMesh(clac3, ptr);

  polyMeshRepartition(&clac3, ptr);

  polyMeshExtractFaces(ptr->GetUnstructCells());

  polyMeshExportImport(&clac3, ptr, MeshPath("output/cube3D_clipped_HexaHexaRotate"), 1, false);
}

TEST(FSClippingTestInterface, ReconstructClippedMesh2Cylinders)
{
  FSClac clac1, clac2, clac3;

  auto mesh1 = LoadMeshTau(clac1, MeshPath("input/mesh_cylinder_1.grid"));
  auto mesh2 = LoadMeshTau(clac2, MeshPath("input/mesh_cylinder_2.grid"));

  const FS_intT marker1 = 1;
  const FS_intT marker2 = 1;
  const FS_floatT tol = 1e-6;
  FSClippingInterface clip(clac1, clac2, tol, marker1, marker2);

  auto meshClipped = clip.BuildVolumeInterface(mesh1, mesh2);
  FSMeshData* ptr = meshClipped.GetMeshData();

  CheckMesh(clac3, ptr);

  polyMeshRepartition(&clac3, ptr);

  polyMeshExtractFaces(ptr->GetUnstructCells());

  polyMeshExportImport(&clac3, ptr, MeshPath("output/2Cylinders"), 1, false);
}

_FS_END_NAMESPACE