#include "FSClipping/FSMeshReconstruction.h"
#include "FSMeshImportParamsTAU.h"
#include "FSMeshExportFilterVTK.h"
#include "FSMeshExportFilterVTK.h"
#include "FSMeshExportParamsVTK.h"
#include "FSMeshCreateLocalNumbering.h"
#include "FSMeshExportFilterHDF5.h"
#include "FSMeshImportFilterHDF5.h"
#include "FSMeshPartitionerRCB.h"
#include "FSMeshPartitionerPARMETIS.h"
#include "gtest/gtest.h"

#include "FSClipping/FSClippingInterface.h"
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSTopologyAssembler.h"
#include "FSClipping/FSClippingFace.h"
#include "FSClipping/FSFaceSeparator.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSClipping/FSPolyFaceBuilder.h"

_FS_BEGIN_NAMESPACE

void polyMeshExportImport(FSClac* clac, FSMeshData* meshDataPtr, const FSString filename_prefix,
                          FS_intT logLevel, bool partitionIndependent)
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

static void
polyMeshExtractFaces(FSUnstructMeshData& unstructMeshData)
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

TEST(FSClippingTestInterface, BuildSurfaceInterface)
{
  FSClac clac1, clac2;
  FSMeshImportParamsTAU params1, params2;

  params1.mMeshFilename =
    //  "/stck/aleprevo/test/test_clipping/hexa_tetra/GRID/cube_hexa.grid";
    //"/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/cube_tetra_coarse.grid";
    "/stck/aleprevo/test/test_clipping/2026-01-27_slidingMeshes-maillage-de-revolution-antonin_mc/interior.grid";
  params2.mMeshFilename =
    //  "/stck/aleprevo/test/test_clipping/hexa_tetra/GRID/cube_hexa_fine.grid";
    //"/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/cube_hexa_coarse.grid";
    "/stck/aleprevo/test/test_clipping/2026-01-27_slidingMeshes-maillage-de-revolution-antonin_mc/exterior_rotate.grid";

  FSMesh mesh1(&clac1);
  FSMesh mesh2(&clac2);
  ASSERT_TRUE(mesh1.ImportMesh(&params1));
  ASSERT_TRUE(mesh2.ImportMesh(&params2));

  mesh1.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  mesh2.GetMeshData()->GetUnstructCells().CreateLocalNumbering();

  FS_intT boundaryMarkerMesh1 = 4;
  FS_intT boundaryMarkerMesh2 = 3;
  FS_floatT tolerance = 1e-10;
  FSClac clac3;

  // compute the interface between mesh1 and mesh2
  FSClippingInterface clippingInterface(clac1, clac2, tolerance, boundaryMarkerMesh1, boundaryMarkerMesh2);
  auto meshDataInterface = clippingInterface.BuildSurfaceInterface(mesh1, mesh2, clac3);
  FSMeshData* meshDataInterfacePtr = &meshDataInterface;

  // print mesh info
  FSMeshPrintInfo printInfo(&clac3);
  FSMeshOpParams dummy;
  bool success = printInfo.DoOp(meshDataInterfacePtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);


  // check mesh
  FSMeshCheck meshCheck(&clac3);
  success = meshCheck.DoOp(meshDataInterfacePtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  FS_intT logLevel = 1;

  polyMeshRepartition(&clac3, meshDataInterfacePtr);
  polyMeshExportImport(&clac3, meshDataInterfacePtr,
                       "/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/Output/rotate_2D", logLevel, false);

  // meshDataInterface.GetUnstructCells().CreateLocalNumbering();
  // polyMeshRepartition(&clac3, meshDataInterfacePtr);
  // polyMeshExportImport(&clac3, meshDataInterfacePtr,
  //                      "/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/Output/cube_hexa_clipped_2D_local", logLevel, false);

  polyMeshExtractFaces(meshDataInterfacePtr->GetUnstructCells());
}


TEST(FSClippingTestInterface, BuildVolumeInterface)
{
  FSClac clac1, clac2;
  FSMeshImportParamsTAU params1, params2;

  params1.mMeshFilename =
    "/stck/aleprevo/test/test_clipping/hexa_tetra/GRID/cube_hexa.grid";
  //"/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/cube_tetra_coarse.grid";
  //"/stck/aleprevo/test/test_clipping/2026-01-27_slidingMeshes-maillage-de-revolution-antonin_mc/interior.grid";
  params2.mMeshFilename =
    "/stck/aleprevo/test/test_clipping/hexa_tetra/GRID/cube_hexa_fine.grid";
  //"/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/cube_hexa_coarse.grid";
  //"/stck/aleprevo/test/test_clipping/2026-01-27_slidingMeshes-maillage-de-revolution-antonin_mc/exterior_rotate.grid";


  FSMesh mesh1(&clac1);
  FSMesh mesh2(&clac2);
  ASSERT_TRUE(mesh1.ImportMesh(&params1));
  ASSERT_TRUE(mesh2.ImportMesh(&params2));

  FSMeshData* meshData1 = mesh1.GetMeshData();
  FSMeshPrintInfo printInfo(&clac1);
  FSMeshOpParams dummy;
  bool success = printInfo.DoOp(meshData1, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  mesh1.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  mesh2.GetMeshData()->GetUnstructCells().CreateLocalNumbering();

  FS_intT boundaryMarkerMesh1 = 2; // 2
  FS_intT boundaryMarkerMesh2 = 1; // 1
  FS_floatT tolerance = 1e-8;
  FSClac clac3;

  // compute the interface between mesh1 and mesh2
  FSClippingInterface clippingInterface(clac1, clac2, tolerance, boundaryMarkerMesh1, boundaryMarkerMesh2);
  auto meshDataInterface = clippingInterface.BuildVolumeInterface(mesh1, mesh2, clac3);
  FSMeshData* meshDataInterfacePtr = &meshDataInterface;

  // print mesh info
  success = printInfo.DoOp(meshDataInterfacePtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);


  // check mesh
  FSMeshCheck meshCheck(&clac3);
  success = meshCheck.DoOp(meshDataInterfacePtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  FS_intT logLevel = 1;

  polyMeshRepartition(&clac3, meshDataInterfacePtr);
  polyMeshExportImport(&clac3, meshDataInterfacePtr,
                       "/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/Output/hexa_clipped", logLevel, false);

  //  meshDataInterface.GetUnstructCells().CreateLocalNumbering();
  //  polyMeshRepartition(&clac3, meshDataInterfacePtr);
  //  polyMeshExportImport(&clac3, meshDataInterfacePtr,
  //                       "/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/Output/hexa_clipped_local", logLevel, false);

  polyMeshExtractFaces(meshDataInterfacePtr->GetUnstructCells());
}
