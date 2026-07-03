#include "FSClipping/FSClippingInterfacePar.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"
#include <FSMeshCheck.h>
#include <FSMeshFaceExtractor.h>
#include <FSMeshPartitionerRCB.h>
#include <FSMeshPrintInfo.h>

_FS_BEGIN_NAMESPACE

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
  FS_intT meshID = globalClac.GetProcID();
  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-8;
  FS_intT marker = -1;

  if(meshID == 1) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine.grid"));
    marker = 1;
  }

  if(meshID == 0) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse.grid"));
    marker = 2;
  }

  FSClippingInterfacePar surfaceInterface(globalClac, clac, tol, marker);
  FSMesh meshClipped = surfaceInterface.BuildSurfaceInterface(mesh);

  if(meshID == 0 || meshID == 1) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    CheckMesh(clac, ptr);
    polyMeshRepartition(&clac, ptr);
  }
}

TEST(FSCLippingTestInterfacePar, VolumeInterface)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT meshID = globalClac.GetProcID();
  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-8;
  FS_intT marker = -1;

  if(meshID == 1) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine.grid"));
    marker = 1;
  }

  if(meshID == 0) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse.grid"));
    marker = 2;
  }

  FSClippingInterfacePar clip(globalClac, clac, tol, marker);
  FSMesh meshClipped = clip.BuildVolumeInterface(mesh);

  if(meshID == 0 || meshID == 1) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    CheckMesh(clac, ptr);
    polyMeshRepartition(&clac, ptr);
  }
}

// The two cylinder boundaries only partially overlap (relative rotation):
// faces at the rim of the interface are not fully covered by the other mesh.
// Both sides must keep those faces whole (NOT_COVERED) so that neither
// rebuilt mesh has holes or missing faces.
TEST(FSCLippingTestInterfacePar, VolumeInterface2Cylinders)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT meshID = globalClac.GetProcID();
  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-6;
  FS_intT marker = -1;

  if(meshID == 0) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/mesh_cylinder_1.grid"));
    marker = 1;
  }

  if(meshID == 1) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/mesh_cylinder_2.grid"));
    marker = 1;
  }

  FSClippingInterfacePar clip(globalClac, clac, tol, marker);
  FSMesh meshClipped = clip.BuildVolumeInterface(mesh);

  if(meshID == 0 || meshID == 1) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    CheckMesh(clac, ptr);
    polyMeshRepartition(&clac, ptr);
  }
}

_FS_END_NAMESPACE
