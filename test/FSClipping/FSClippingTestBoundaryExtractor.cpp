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

TEST(FSClippingTestBoundaryFaceProvider, ComputeGeometry)
{
  FSClac clac;
  auto mesh = LoadMesh(clac, MeshPath("input/cube_hexa_coarse.grid"));

  FSMeshFaceExtractor ex;
  auto boundaryExtraction = FSBoundaryFaceProvider::Extract(mesh, ex, 1, 1e-12);
  auto clipFaces = boundaryExtraction.faces;
  for(const auto& f : clipFaces) {
    // tests géométriques
    EXPECT_NEAR(f.tangentVector1().L2Norm(), 1.0, 1e-12);
    EXPECT_NEAR(f.normal().L2Norm(), 1.0, 1e-12);
  }
}

TEST(FSClippingTestBoundaryFaceProvider, GeomFaceKeyCubeHexaCoarse)
{
  FSClac clac;
  auto mesh = LoadMesh(clac, MeshPath("input/cube_hexa_coarse.grid"));

  FSMeshFaceExtractor ex;
  auto boundaryExtraction = FSBoundaryFaceProvider::Extract(mesh, ex, 1, 1e-6);

  EXPECT_EQ(boundaryExtraction.faces.size(), 4);

  // All faces should be stored in the faceKeys set
  EXPECT_EQ(boundaryExtraction.faceKeys.size(), 4);

  // Verify each face key has 4 vertices (quad faces)
  for(const auto& key : boundaryExtraction.faceKeys) {
    EXPECT_EQ(key.pts.size(), 4);
  }

  // Verify that inserting the same key again doesn't increase the set size
  const auto& firstKey = *boundaryExtraction.faceKeys.begin();
  auto result = boundaryExtraction.faceKeys.insert(firstKey);
  EXPECT_FALSE(result.second); // insertion should fail (already exists)
  EXPECT_EQ(boundaryExtraction.faceKeys.size(), 4);
}

_FS_END_NAMESPACE
