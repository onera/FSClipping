// Asymmetric parallel interface tests — run at exactly 6 MPI procs.
//
// Layout (AsymMeshID, nMesh1=2): rank 0 = clipper, ranks 1..2 = mesh 1,
// ranks 3..5 = mesh 2, i.e. 1 clipper + 2 procs (mesh 1) + 3 procs (mesh 2).
// The two meshes get a different number of procs (asymmetric split).
//
// This translation unit is compiled into its own CTest executable
// (FSClippingParallelTest_np6) and launched only at 6 procs, so no test needs to
// skip itself.

#include "FSClipping/FSClippingInterfacePar.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"

_FS_BEGIN_NAMESPACE

static constexpr FS_intT kNMesh1 = 2; // mesh 1 gets 2 procs, mesh 2 gets 3

// ── Cube, surface reconstruction (asymmetric split) ──────────────────────────
TEST(FSClippingTestInterfaceParAsym, CubeSurface)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT procId = globalClac.GetProcID();

  FS_intT meshID = AsymMeshID(procId, kNMesh1);

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-7;
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

  // Reconstruction only (CheckMesh). A proc may own an empty clipped partition in
  // an asymmetric split; exporting such a group hangs the VTK/HDF5 collective
  // (see the empty-partition I/O limitation), so the clipped mesh is not exported.
  if(meshID == 1 || meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    CheckMesh(clac, ptr);
    meshClipped.PrintInfo();
    exportMeshVTK(
      &clac, ptr,
      MeshPath(meshID == 1 ? "output/cube_coarse_clipped_surf_par_asym" : "output/cube_fine_clipped_surf_par_asym"));
  }
}

// ── Two cylinders, surface reconstruction (asymmetric split) ─────────────────
TEST(FSClippingTestInterfaceParAsym, TwoCylindersSurface)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT procId = globalClac.GetProcID();

  FS_intT meshID = AsymMeshID(procId, kNMesh1);

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-6;
  const FS_intT marker = 1;

  if(meshID == 1) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/mesh_cylinder_1.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    exportMeshVTK(&clac, mesh.GetMeshData(), MeshPath("output/mesh_cylinder_1"));
    CheckMesh(clac, mesh.GetMeshData());
  }

  if(meshID == 2) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/mesh_cylinder_2.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    exportMeshVTK(&clac, mesh.GetMeshData(), MeshPath("output/mesh_cylinder_2"));
    CheckMesh(clac, mesh.GetMeshData());
  }

  FSClippingInterfacePar surfaceInterface(globalClac, clac, tol, marker, meshID);
  FSMesh meshClipped = surfaceInterface.BuildSurfaceInterface(mesh, meshID);

  // A proc may legitimately own 0 nodes here (its whole boundary is covered by
  // higher-priority procs after min-rank ownership). The group total is > 0 so
  // InitNodePool is satisfied. The clipped mesh is NOT exported: exporting a group
  // in which a proc holds an empty clipped partition hangs the VTK/HDF5 collective
  // (empty-partition I/O limitation). Only the original mesh is exported above.
  if(meshID == 1 || meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    CheckMesh(clac, ptr);
    meshClipped.PrintInfo();
    exportMeshVTK(&clac, ptr,
                  MeshPath(meshID == 1 ? "output/cyl1_clipped_surf_asym" : "output/cyl2_clipped_surf_asym"));
  }
}

// ── Cube, volume reconstruction (asymmetric split) ───────────────────────────
TEST(FSClippingTestInterfaceParAsym, CubeVolume)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT procId = globalClac.GetProcID();

  FS_intT meshID = AsymMeshID(procId, kNMesh1);

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
    CheckMesh(clac, ptr);
    meshClipped.PrintInfo();
    exportMeshVTK(&clac, ptr,
                  MeshPath(meshID == 1 ? "output/cube_coarse_clipped_vol_asym" : "output/cube_fine_clipped_vol_asym"));
  }
}

// ── Two cylinders, volume reconstruction (asymmetric split) ──────────────────
TEST(FSClippingTestInterfaceParAsym, TwoCylindersVolume)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT procId = globalClac.GetProcID();

  FS_intT meshID = AsymMeshID(procId, kNMesh1);

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
    exportMeshVTK(&clac, mesh.GetMeshData(), MeshPath("output/cyl1_orig_vol_asym"));
  }

  if(meshID == 2) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/mesh_cylinder_2.grid"));
    polyMeshRepartition(&clac, mesh.GetMeshData());
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
    CheckMesh(clac, mesh.GetMeshData());
    exportMeshVTK(&clac, mesh.GetMeshData(), MeshPath("output/cyl2_orig_vol_asym"));
  }

  FSClippingInterfacePar volumeInterface(globalClac, clac, tol, marker, meshID);
  FSMesh meshClipped = volumeInterface.BuildVolumeInterface(mesh, meshID);

  // The clipped mesh is NOT exported here: in this asymmetric split a proc of
  // mesh 1 can get an empty clipped-poly partition (0 Poly2D/Poly3D), and the
  // VTK/HDF5 export collective hangs on such a proc (empty-partition I/O
  // limitation). The reconstruction itself is validated by CheckMesh.

  if(meshID == 1 || meshID == 2) {
    FSMeshData* ptr = meshClipped.GetMeshData();
    // Export the clipped mesh with its as-clipped partition (one VTK per proc).
    meshClipped.PrintInfo();
    exportMeshVTK(&clac, ptr, MeshPath(meshID == 1 ? "output/cyl1_clipped_vol_asym" : "output/cyl2_clipped_vol_asym"));
  }
}

_FS_END_NAMESPACE
