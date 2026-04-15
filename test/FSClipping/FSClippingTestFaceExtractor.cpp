#include "FSMeshImportParamsTAU.h"

#include "gtest/gtest.h"

#include "FSClipping/FSClippingFace.h"
#include "FSClipping/FSFaceSeparator.h"

#include "TestUtils.hpp"


_FS_BEGIN_NAMESPACE


FSMesh LoadMesh(FSClac& clac, const std::string& filename)
{
  FSMeshImportParamsTAU params;
  params.mMeshFilename = filename;

  FSMesh mesh(&clac);
  EXPECT_TRUE(mesh.ImportMesh(&params));

  mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  return mesh;
}

TEST(FSClippingTestFaceExtractor, ComputeGeometry)
{
  FSClac clac;
  auto mesh = LoadMesh(clac, MeshPath("input/cube_hexa_coarse.grid"));

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
