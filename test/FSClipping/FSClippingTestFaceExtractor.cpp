#include "FSMeshImportParamsTAU.h"

#include "gtest/gtest.h"

#include "FSClipping/FSClippingFace.h"
#include "FSClipping/FSFaceSeparator.h"

_FS_BEGIN_NAMESPACE

TEST(FSClippingTestFaceExtractor, ComputeGeometry)
{
  FSClac clac;
  FSMeshImportParamsTAU params;
  // params.mMeshFilename = "${HOME}/path/to/hexa.grid";
  params.mMeshFilename =
    "${HOME}/code_dev/CODA_src/FSClipping/test/Mesh/hexa.grid";
  // hexa.grid is in FSClipping/test/Mesh
  // but you can load any mesh.grid mesh file
  FSMesh mesh(&clac);
  ASSERT_TRUE(mesh.ImportMesh(&params));

  FSMeshFaceExtractor ex;
  mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  ASSERT_TRUE(ex.PrepareFaceConnectivity(mesh.GetMeshData()->GetUnstructCells(),
                                         true, false));
  ASSERT_TRUE(ex.PrepareFaceNodeCoordinates(FSQuantityDescArrayT()));

  auto bdryFaces =
    FSFaceSeparator::SeparateBoundariesFaceWithMarker(mesh, ex, 1);
  ASSERT_FALSE(bdryFaces.empty());

  std::vector<FSClippingFace> clipFaces;
  FSFloatArrayT faceNodeCoordinates;
  for(const auto& f : bdryFaces) {
    const FS_intT faceIndex = ex.GetFaceIndex(*(f._faceFSDM));
    ex.GetFaceNodeCoordinates(faceIndex, faceNodeCoordinates);
    clipFaces.emplace_back(f, faceIndex, faceNodeCoordinates);

    // tests géométriques
    EXPECT_NEAR(clipFaces.back().tangentVector1().L2Norm(), 1.0, 1e-12);
    EXPECT_NEAR(clipFaces.back().normal().L2Norm(), 1.0, 1e-12);
  }
}

// TEST(FSClippingTestFaceExtractor, ExtractAttribute)
//{
//   FSClac clac;
//   FSMeshImportParamsTAU params;
//   // params.mMeshFilename = "${HOME}/path/to/hexa.grid";
//   params.mMeshFilename =
//     "${HOME}/code_dev/CODA_src/FSClipping/test/Mesh/hexa.grid";
//   // hexa.grid is in FSClipping/test/Mesh
//   // but you can load any mesh.grid mesh file
//   FSMesh mesh(&clac);
//   ASSERT_TRUE(mesh.ImportMesh(&params));
//
//   FSMeshData* meshData = mesh.GetMeshData();
//   FSUnstructMeshData& unstructCell = meshData->GetUnstructCells();
//   FSMeshPrintInfo printInfo(&clac);
//   FSMeshOpParams dummy;
//   bool success = printInfo.DoOp(meshData, &dummy);
//   if(!(success)) {
//     FSError.Print();
//   }
//   ASSERT_TRUE(success);
//
//   const FSString attrCADGroupIDName =
//     FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
//   auto cellTypes = mesh.GetCellTypes();
//   for(const auto& t : cellTypes) {
//     const auto& cellAttribNames = unstructCell.GetCellAttributes(t);
//
//     for(FSStringArrayT::ConstIterator AI = cellAttribNames.BeginConst(); AI.IsValid(); ++AI) {
//       const FSIntArrayT& valuesOrig = unstructCell.GetCellAttribute(*AI, t);
//       std::cout << *AI << " : " << valuesOrig.Size(0) << std::endl;
//       for(FS_intT i = 0; i < valuesOrig.Size(); i++) {
//         std::cout << valuesOrig[i] << " " << "\t";
//       }
//       std::cout << "\n";
//     }
//   }
//  for(FSStringArrayT::ConstIterator AI = )

/* for(const auto& t : cellTypes) {
   std::cout << FSMeshEnums::CellTypeToString(t) << " : ";
   if(FSMeshEnums::IsUnstructSurfaceCellType(t)) {
     auto cellPool = unstructCell.GetCellPool(t);

     if(cellPool->HasCellAttribute(attrCADGroupIDName)) {
       const auto& boundaryMarkersValue = cellPool->GetCellAttribute(attrCADGroupIDName);
       std::cout << boundaryMarkersValue.Size() << std::endl;
       std::cout << cellPool->GetNAttributes() << std::endl;
       std::cout << cellPool->GetNCellAttributes(attrCADGroupIDName) << std::endl;

       for(FS_intT i = 0; i < boundaryMarkersValue.Size(); i++) {
         std::cout << boundaryMarkersValue[i] << "\t";
       }
     }
   }
   std::cout << std::endl;
 }
}*/
_FS_END_NAMESPACE
