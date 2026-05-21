#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSClippingInterfacePar.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSTopologyAssembler.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"

_FS_BEGIN_NAMESPACE

TEST(FSCLippingTestInterfaceParllel3, SurfaceInterface)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT meshID = globalClac.GetProcID();
  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-8;
  FS_intT marker = -1;

  if(meshID == 0) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse.grid"));
    marker = 2;

  } else if(meshID == 1) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine.grid"));
    marker = 1;
  }

  FSClippingInterfacePar surfaceInterface(globalClac, clac, tol, marker);
  surfaceInterface.BuildSurfaceInterface(mesh);
}

_FS_END_NAMESPACE
