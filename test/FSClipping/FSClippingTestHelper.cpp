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
  // FSMeshPrintInfo printInfo(clac);
  FSMeshOpParams dummy;
  // success = printInfo.DoOp(meshDataImpPtr, &dummy);
  // if(!(success)) {
  //   FSError.Print();
  // }
  // ASSERT_TRUE(success);

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
  // FSMeshPrintInfo printInfo(clac);
  //  repartition of mesh
  FSMeshPartitionerRCB partitioner(clac);
  FSMeshPartitioningParamsRCB partitionParams;
  success = partitioner.DoOp(meshDataPtr, &partitionParams);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  // print mesh info
  // success = printInfo.DoOp(meshDataPtr, &dummy);
  // if(!(success)) {
  //  FSError.Print();
  //}
  // ASSERT_TRUE(success);

  // check mesh
  FSMeshCheck meshCheck(clac);
  success = meshCheck.DoOp(meshDataPtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);
}

TEST(FSClippingTestHelper, Interface)
{
  FSClac clac1, clac2;
  FSMeshImportParamsTAU params1, params2;
  // params.mMeshFilename = "${HOME}/path/to/hexa.grid";
  params1.mMeshFilename =
    "/stck/aleprevo/test/test_clipping/2026-01-19_slidingMeshes-boites-antonin_mc/cube_tetra.grid";
  //"/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/cube_hexa_coarse.grid";
  // params.mMeshFilename = "${HOME}/path/to/hexa.grid";
  params2.mMeshFilename =
    "/stck/aleprevo/test/test_clipping/2026-01-19_slidingMeshes-boites-antonin_mc/cube_quad.grid";
  //"/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/cube_tetra_coarse.grid";

  // hexa.grid is in FSClipping/test/Mesh
  // but you can load any mesh.grid mesh file
  FSMesh mesh1(&clac1);
  FSMesh mesh2(&clac2);
  ASSERT_TRUE(mesh1.ImportMesh(&params1));
  ASSERT_TRUE(mesh2.ImportMesh(&params2));

  FSMeshFaceExtractor ex1, ex2;
  mesh1.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  ex1.PrepareFaceConnectivity(mesh1.GetMeshData()->GetUnstructCells(), true,
                              false);
  ex1.PrepareFaceNodeCoordinates(FSQuantityDescArrayT());

  mesh2.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  ex2.PrepareFaceConnectivity(mesh2.GetMeshData()->GetUnstructCells(), true,
                              false);
  ex2.PrepareFaceNodeCoordinates(FSQuantityDescArrayT());

  auto bdryFaces1 =
    FSFaceSeparator::SeparateBoundariesFaceWithMarker(mesh1, ex1, 1);

  auto bdryFaces2 =
    FSFaceSeparator::SeparateBoundariesFaceWithMarker(mesh2, ex2, 2);

  FSQuantityDescArrayT coordDesc;
  FS_intT nodeDatasetOffset = 0;
  FSFloatArrayT coords;

  mesh1.GetMeshData()->GetUnstructCells().GetCoordinates3D(coordDesc, coords,
                                                           nodeDatasetOffset);
  std::vector<FSClippingFace> clipFaces1, clipFaces2;

  FSFloatArrayT faceCoordinates1;
  FSFloatArrayT faceCoordinates2;
  for(const auto& f : bdryFaces1) {
    FS_intT faceIndex = ex1.GetFaceIndex(*(f._faceFSDM));
    ex1.GetFaceNodeCoordinates(faceIndex, faceCoordinates1);
    clipFaces1.emplace_back(f, faceIndex, faceCoordinates1);
  }
  for(const auto& f : bdryFaces2) {
    FS_intT faceIndex = ex2.GetFaceIndex(*(f._faceFSDM));
    ex2.GetFaceNodeCoordinates(faceIndex, faceCoordinates2);
    clipFaces2.emplace_back(f, faceIndex, faceCoordinates2);
  }
  FSFaceMatcher matcher(clac2, clipFaces1, clipFaces2, 1e-6);

  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);
  // ASSERT_TRUE(matches.size() == 16);

  // std::unordered_map<FS_intT, std::vector<FSVec3> > cell2NodeNew;
  FSCell2NodeBuilder cell2NodeBuilder(1e-6);
  // std::vector<FSVec3> newCord;


  for(const auto& f : matches) {
    // retrouver bdr face associe
    auto elemIndex = f.elemOwner1; //-> indice dans FSDM
    auto cellType = f.elemOwnerType1;
    cell2NodeBuilder.AddOldCellNodes(elemIndex, cellType);

    if(f.type != FSFaceMatch::UNKNOWN) {
      cell2NodeBuilder.AddClippedPolygon(elemIndex, f.clippedPoly3D);
    }
  }

  cell2NodeBuilder.BuildGlobalNumbering();

  auto globalCoords = cell2NodeBuilder.GlobalCoords();
  auto cellData = cell2NodeBuilder.CellData();
  auto cell2Node = cell2NodeBuilder.Cell2Node();

  // for(const auto& p : cellData) {
  //   std::cout << p.first << " : ";
  //   for(const auto& c : p.second.coords)
  //     std::cout << "(" << c[0] << ", " << c[1] << ", " << c[2] << ")" << "; ";
  //   std::cout << "\n";
  // }



  // for(const auto& p : cellData) {
  //   std::cout << p.first << " : ";
  //   for(const auto& n : p.second.nodeIds)
  //     std::cout << n << "; ";
  //   std::cout << "\n";
  // }


  // std::cout << "\n"
  //           << "Global coords : ";
  // for(const auto& c : globalCoords)
  //   std::cout << "(" << c[0] << ", " << c[1] << ", " << c[2] << ")" << " ";

  // std::cout << "\n";


  FSMeshPolyFaceStorage polyFaces;
  FSPolyFaceBuilder polyFaceBuilder(matches, cell2NodeBuilder, 1e-6);
  polyFaceBuilder.Build(polyFaces);
  //  auto cell2Node = cell2NodeBuilder.cell2Node();

  FSClac clac3;
  FSMesh meshInterface = FSMesh(&clac3);
  FSMeshData meshDataInterface = FSMeshData(&clac3);
  // init pointer to meshdata
  FSMeshData* meshDataInterfacePtr = &meshDataInterface;
  // init reference to unstructured mesh data
  FSUnstructMeshData& unstructMeshData = meshDataInterface.GetUnstructCells();


  // // begin initialization
  meshDataInterface.BeginInitialization();

  meshDataInterface.InitUnstructNodes(cell2NodeBuilder.GlobalCoords().size());
  meshDataInterface.InitUnstructCells(FSMeshEnums::CT_Poly3D, cell2Node);
  meshDataInterface.InitUnstructCellFaces(FSMeshEnums::CT_Poly3D, polyFaces);


  FSQuantityDescArrayT coordsDesc(3);
  coordsDesc[0] = FSQuantityDesc(FSDataName::Coordinates(), FSDataName::Coordinate().X());
  coordsDesc[1] = FSQuantityDesc(FSDataName::Coordinates(), FSDataName::Coordinate().Y());
  coordsDesc[2] = FSQuantityDesc(FSDataName::Coordinates(), FSDataName::Coordinate().Z());

  FSFloatArrayT coordInterface(globalCoords.size(), 3);
  for(std::size_t i = 0; i < globalCoords.size(); i++) {
    coordInterface(i, 0) = globalCoords[i][0];
    coordInterface(i, 1) = globalCoords[i][1];
    coordInterface(i, 2) = globalCoords[i][2];
  }

  bool success = unstructMeshData.SetCoordinates3D(coordsDesc, coordInterface);
  ASSERT_TRUE(success);


  FS_intT currentOffset = 0;
  FSIntArrayT cellTypeArray_2 = meshDataInterface.GetUnstructCells().GetCellTypesArray();
  for(FSIntArrayT::ConstIterator cellType = cellTypeArray_2.BeginConst(); cellType.IsValid(); cellType.Next()) {
    meshDataInterface.GetUnstructCells().InitGlobalCellNumber((FSMeshEnums::CellType)*cellType, currentOffset);
    currentOffset += meshDataInterface.GetUnstructCells().GetNCells((FSMeshEnums::CellType)*cellType);
  }

  // end initialization
  meshDataInterface.EndInitialization();
  // check if initialization is complete
  success = meshDataInterface.IsInitialized();
  ASSERT_TRUE(success);

  // print mesh info
  // FSMeshPrintInfo printInfo(&clac3);
  FSMeshOpParams dummy;
  // success = printInfo.DoOp(meshDataInterfacePtr, &dummy);
  // if(!(success)) {
  //   FSError.Print();
  // }
  // ASSERT_TRUE(success);


  // // // check mesh
  FSMeshCheck meshCheck(&clac3);
  success = meshCheck.DoOp(meshDataInterfacePtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  FS_intT logLevel = 1;

  polyMeshRepartition(&clac3, meshDataInterfacePtr);
  polyMeshExportImport(&clac3, meshDataInterfacePtr, "/stck/aleprevo/FSTwoPolys2DGlobal", logLevel, false);
  polyMeshExportImport(&clac3, meshDataInterfacePtr, "FSTwoPolys2DGlobal_partind", logLevel, true);


  meshDataInterface.GetUnstructCells().CreateLocalNumbering();
  polyMeshRepartition(&clac3, meshDataInterfacePtr);
  polyMeshExportImport(&clac3, meshDataInterfacePtr, "/stck/aleprevo/FSTwoPolys2DLocal", logLevel, false);
  polyMeshExportImport(&clac3, meshDataInterfacePtr, "FSTwoPolys2DLocal_partind", logLevel, true);
}
_FS_END_NAMESPACE
