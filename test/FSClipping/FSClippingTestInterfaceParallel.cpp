#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSClippingInterface.h"
#include "TestUtils.hpp"
#include "gtest/gtest.h"

#include "FSMeshImportParamsTAU.h"
#include <FSMesh.h>

_FS_BEGIN_NAMESPACE

namespace {

FSMesh LoadMeshWithClac(FSClac& clac, const std::string& filename)
{
  FSMeshImportParamsTAU params;
  params.mMeshFilename = filename;

  FSMesh mesh(&clac);
  EXPECT_TRUE(mesh.ImportMesh(&params));
  mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  return mesh;
}

} // namespace

// Two-process test: each proc owns one mesh and extracts its boundary faces independently.
// Proc 0 → cube_hexa_coarse (marker 2), Proc 1 → cube_hexa_fine (marker 1).
TEST(FSClippingTestInterfaceParallel, ExtractFacesPerProc)
{
  // Split MPI_COMM_WORLD so each proc gets its own 1-proc communicator.
  FSClac worldClac(MPI_COMM_WORLD);
  FSClac localClac;
  worldClac.DivideIntoGroups(worldClac.GetProcID(), localClac);

  const FS_floatT tol = 1e-8;

  BoundaryExtraction extraction;

  if(worldClac.GetProcID() == 0) {
    auto mesh = LoadMeshWithClac(localClac, MeshPath("input/cube_hexa_coarse.grid"));
    const FS_intT marker = 2;
    FSClippingInterface clip(localClac, localClac, tol, marker, marker);
    extraction = clip.BuildSurfaceInterfacePar(mesh, marker);
  } else {
    auto mesh = LoadMeshWithClac(localClac, MeshPath("input/cube_hexa_fine.grid"));
    const FS_intT marker = 1;
    FSClippingInterface clip(localClac, localClac, tol, marker, marker);
    extraction = clip.BuildSurfaceInterfacePar(mesh, marker);
  }

  EXPECT_GT(extraction.faces.size(), 0u);
}

_FS_END_NAMESPACE
