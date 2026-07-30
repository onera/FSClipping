#pragma once

#include "../TestUtils.hpp"
#include "FSMeshImportParamsTAU.h"
#include "gtest/gtest.h"
#include <FSClac.h>
#include <FSMesh.h>
#include <FSMeshCheck.h>
#include <FSMeshData.h>
#include <FSMeshFaceExtractor.h>
#include <FSMeshExportFilterHDF5.h>
#include <FSMeshExportFilterVTK.h>
#include <FSMeshExportParamsHDF5.h>
#include <FSMeshImportFilterHDF5.h>
#include <FSMeshPartitionerRCB.h>
#include <FSMeshPrintInfo.h>

_FS_BEGIN_NAMESPACE

inline FSMesh LoadMeshWithClac(FSClac& clac, const FSString& filename)
{
  FSMeshImportParamsTAU params;
  params.mMeshFilename = filename;
  FSMesh mesh(&clac);
  EXPECT_TRUE(mesh.ImportMesh(&params));
  return mesh;
}

// ── Proc-layout helpers ──────────────────────────────────────────────────────
//
// Every parallel interface test is now compiled into exactly one CTest
// executable, run at exactly the MPI world size its layout needs (np5 for the
// symmetric 1+2+2 layout, np6 for the asymmetric 1+2+3 layout). There is no
// cross-proc-count registration any more, so tests no longer skip themselves —
// the SKIP_UNLESS_* macro below is kept only for ad-hoc / manual runs.

#define SKIP_UNLESS_EXACT_PROCS(clac, n)                                                                               \
  do {                                                                                                                 \
    const FS_intT nProcs__ = (clac).GetNProcs();                                                                       \
    if(nProcs__ != (n))                                                                                                \
      GTEST_SKIP() << "needs exactly " << (n) << " MPI procs, running with " << nProcs__;                              \
  } while(0)

// Symmetric parity layout: rank 0 = clipper, even ranks = mesh 1, odd ranks =
// mesh 2. With np5 both meshes get 2 procs.
inline FS_intT SymMeshID(FS_intT procId)
{
  if(procId == 0)
    return 0; // Clipper
  return (procId % 2 == 0) ? 1 : 2;
}

// Asymmetric layout: rank 0 = clipper, the next nMesh1 ranks build mesh 1, the
// remaining ranks build mesh 2. Contiguous ranges (not rank parity) so the two
// meshes can have a different number of procs. With np6 and nMesh1=2 the split
// is 1 clipper + 2 (mesh 1) + 3 (mesh 2). Returns meshID in {0,1,2}.
inline FS_intT AsymMeshID(FS_intT procId, FS_intT nMesh1)
{
  if(procId == 0)
    return 0; // Clipper
  if(procId <= nMesh1)
    return 1;
  return 2;
}

// ── Mesh-op helpers shared by the interface tests ────────────────────────────

// Export the mesh as one VTK (legacy ASCII) file per proc: prefix_<procID>.vtk.
// Used to visualise the per-proc partition, both the original mesh (before
// clipping) and the reconstructed mesh (after clipping). No HDF5, no import.
inline void exportMeshVTK(FSClac* clac, FSMeshData* meshDataPtr, const FSString& filenamePrefix)
{
  const FS_intT procID = FSCLAC_PROCID(clac);
  FSMeshExportFilterVTK exportFilterVTK(clac);
  FSMeshExportParamsVTK exportParamsVTK;
  exportParamsVTK.mMeshFilename = filenamePrefix + FSString(".vtk");
  exportParamsVTK.mFormatType = FSVtkEnums::FT_Raw;
  exportParamsVTK.mFilePerProcess = false;
  bool success = exportFilterVTK.DoOp(meshDataPtr, &exportParamsVTK);
  if(!success)
    FSError.Print();
  ASSERT_TRUE(success);
}

inline void polyMeshRepartition(FSClac* clac, FSMeshData* meshDataPtr)
{
  FSMeshOpParams dummy;
  bool success;
  //  repartition of mesh
  FSMeshPartitionerRCB partitioner(clac);
  FSMeshPartitioningParamsRCB partitionParams;
  success = partitioner.DoOp(meshDataPtr, &partitionParams);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);
}

// Extract the boundary faces and their node coordinates from the mesh. Same
// helper as in the sequential interface test; run before polyMeshExportImport.
inline void polyMeshExtractFaces(FSUnstructMeshData& unstructMeshData)
{
  bool success;
  FSMeshFaceExtractor fex;
  // don't match remote faces, keep pseudo-cells
  success = fex.PrepareFaceConnectivity(unstructMeshData, false, false);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  success = fex.PrepareFaceNodeCoordinates(FSQuantityDescArrayT());
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  FS_intT nFaces = fex.GetNLocalFaces();
  FSIntArrayT faces(nFaces, 6);
  success = fex.ExtractFaces(unstructMeshData, faces, false);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  for(FS_intT face = 0; face < nFaces; ++face) {
    const FS_intT maxFaceNodes = 7; // Maximum number of nodes that we expect of a face
    FSFloatArrayT outFaceNodeCoordinates(maxFaceNodes, 3);
    (void)fex.GetFaceNodeCoordinates(face, outFaceNodeCoordinates);
  }
}

inline void CheckMesh(FSClac& clac, FSMeshData* meshData)
{
  FSMeshOpParams dummy;
  FSMeshCheck meshCheck(&clac);

  ASSERT_TRUE(meshCheck.DoOp(meshData, &dummy)) << "Check failed";
}

inline void polyMeshExportImport(FSClac* clac, FSMeshData* meshDataPtr, const FSString filename_prefix,
                                 FS_intT logLevel, bool partitionIndependent)
{
  FS_intT nProcs = FSCLAC_NPROCS(clac);
  FS_intT procID = FSCLAC_PROCID(clac);
  bool success;

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
  FSMeshData meshDataImp = FSMeshData(clac);
  FSMeshData* meshDataImpPtr = &meshDataImp;
  success = importFilterHDF5.DoOp(meshDataImpPtr, &importParamsHDF5);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);

  FSLog(clac, procID, "################## Check Mesh ###################################\n", logLevel);
  FSMeshCheck meshCheck(clac);
  FSMeshOpParams dummy;
  success = meshCheck.DoOp(meshDataImpPtr, &dummy);
  if(!(success)) {
    FSError.Print();
  }
  ASSERT_TRUE(success);
  // --- end import mesh ---
}

_FS_END_NAMESPACE
