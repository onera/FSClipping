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

_FS_END_NAMESPACE
