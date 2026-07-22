// Symmetric parallel interface tests — run at exactly 5 MPI procs.
//
// Layout (SymMeshID): rank 0 = clipper, even ranks = mesh 1, odd ranks = mesh 2,
// i.e. 1 clipper + 2 procs (mesh 1) + 2 procs (mesh 2). Both meshes get the same
// number of procs.
//
// This translation unit is compiled into its own CTest executable
// (FSClippingParallelTest_np5) and launched only at 5 procs, so no test needs to
// skip itself.

#include "FSClipping/FSClippingInterfacePar.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"

_FS_BEGIN_NAMESPACE

// ── Cube, surface reconstruction ─────────────────────────────────────────────
TEST(FSClippingTestInterfaceParSym, CubeSurface)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT procId = globalClac.GetProcID();

  FS_intT meshID = SymMeshID(procId);

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-8;
  FS_intT marker = -1;

  if(meshID == 1) {
    marker = 6;
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse_par.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, mesh.GetMeshData());
  }

  if(meshID == 2) {
    marker = 5;
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine_par.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, mesh.GetMeshData());
  }

  FSClippingInterfacePar surfaceInterface(globalClac, clac, tol, marker, meshID);
  FSMesh meshClipped = surfaceInterface.BuildSurfaceInterface(mesh, meshID);

  if(meshID == 1 || meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshRepartition(&clac, ptr);
    ptr->GetUnstructCells().CreateLocalNumbering();
    exportMeshVTK(&clac, ptr,
                  MeshPath(meshID == 1 ? "output/cube_coarse_clipped_surf_par" : "output/cube_fine_clipped_surf_par"));
  }
  if(meshID == 1) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshExtractFaces(ptr->GetUnstructCells());
    polyMeshExportImport(&clac, ptr, MeshPath("output/cube_clipped_coarse_par"), 1, false);
  }
  if(meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshExtractFaces(ptr->GetUnstructCells());
    polyMeshExportImport(&clac, ptr, MeshPath("output/cube_clipped_fine_par"), 1, false);
  }
}

// ── Cube, volume reconstruction ──────────────────────────────────────────────
TEST(FSClippingTestInterfaceParSym, CubeVolume)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT procId = globalClac.GetProcID();

  FS_intT meshID = SymMeshID(procId);

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-8;
  FS_intT marker = -1;

  if(meshID == 1) {
    marker = 6;
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse_par.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, mesh.GetMeshData());
  }

  if(meshID == 2) {
    marker = 5;
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine_par.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, mesh.GetMeshData());
  }

  FSClippingInterfacePar volumeInterface(globalClac, clac, tol, marker, meshID);
  FSMesh meshClipped = volumeInterface.BuildVolumeInterface(mesh, meshID);

  if(meshID == 1 || meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshRepartition(&clac, ptr);
    ptr->GetUnstructCells().CreateLocalNumbering();
    exportMeshVTK(&clac, ptr,
                  MeshPath(meshID == 1 ? "output/cube_coarse_clipped_vol_par" : "output/cube_fine_clipped_vol_par"));
  }
  if(meshID == 1) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshExtractFaces(ptr->GetUnstructCells());
    polyMeshExportImport(&clac, ptr, MeshPath("output/cube_clipped_vol_coarse_par"), 1, true);
  }
  if(meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshExtractFaces(ptr->GetUnstructCells());
    polyMeshExportImport(&clac, ptr, MeshPath("output/cube_clipped_vol_fine_par"), 1, true);
  }
}

// ── Two cylinders, volume reconstruction ─────────────────────────────────────
// The two cylinder boundaries only partially overlap (relative rotation): faces
// at the rim of the interface are not fully covered by the other mesh. Both sides
// must keep those faces whole (NOT_COVERED) so that neither rebuilt mesh has
// holes or missing faces.
TEST(FSClippingTestInterfaceParSym, TwoCylindersVolume)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT procId = globalClac.GetProcID();

  FS_intT meshID = SymMeshID(procId);

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-6;
  const FS_intT marker = 1;

  if(meshID == 1) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/mesh_cylinder_1.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, mesh.GetMeshData());
  }

  if(meshID == 2) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/mesh_cylinder_2.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, mesh.GetMeshData());
  }

  FSClippingInterfacePar volumeInterface(globalClac, clac, tol, marker, meshID);
  FSMesh meshClipped = volumeInterface.BuildVolumeInterface(mesh, meshID);

  if(meshID == 1 || meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    // Export the clipped mesh with its as-clipped partition (one VTK per proc).
    polyMeshRepartition(&clac, ptr);
    ptr->GetUnstructCells().CreateLocalNumbering();
    exportMeshVTK(&clac, ptr, MeshPath(meshID == 1 ? "output/cyl1_clipped_vol_view" : "output/cyl2_clipped_vol_view"));
  }
  if(meshID == 1) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshExtractFaces(ptr->GetUnstructCells());
    polyMeshExportImport(&clac, ptr, MeshPath("output/cyl1_clipped_vol_par"), 1, true);
  }
  if(meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    polyMeshExtractFaces(ptr->GetUnstructCells());
    polyMeshExportImport(&clac, ptr, MeshPath("output/cyl2_clipped_vol_par"), 1, true);
  }
}

_FS_END_NAMESPACE
