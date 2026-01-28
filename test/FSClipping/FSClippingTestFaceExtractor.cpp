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

_FS_END_NAMESPACE
