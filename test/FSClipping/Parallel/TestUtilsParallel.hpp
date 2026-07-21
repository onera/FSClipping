#pragma once

#include "../TestUtils.hpp"
#include "FSMeshImportParamsTAU.h"
#include "gtest/gtest.h"
#include <FSClac.h>
#include <FSMesh.h>
#include <FSMeshData.h>

_FS_BEGIN_NAMESPACE

inline FSMesh LoadMeshWithClac(FSClac& clac, const FSString& filename)
{
  FSMeshImportParamsTAU params;
  params.mMeshFilename = filename;
  FSMesh mesh(&clac);
  EXPECT_TRUE(mesh.ImportMesh(&params));
  // mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  return mesh;
}

// A given parallel test is written for a specific MPI world size (its proc
// layout). The same binary is registered in CTest at several proc counts (3, 5,
// 7), so each test must bow out cleanly when the current size is not the one it
// needs — otherwise it would deadlock or segfault. Call these at the very top of
// a TEST body, before any collective. Note: SKIP still lets the test exit its
// body, so no collective must run before the skip check.

// Skip unless the world size is exactly n (fixed-layout tests, e.g. the 3-proc
// send/receive tests, or the 7-proc asymmetric 4+2+1 layout).
#define SKIP_UNLESS_EXACT_PROCS(clac, n)                                                                                \
  do {                                                                                                                  \
    const FS_intT nProcs__ = (clac).GetNProcs();                                                                        \
    if(nProcs__ != (n))                                                                                                 \
      GTEST_SKIP() << "needs exactly " << (n) << " MPI procs, running with " << nProcs__;                               \
  } while(0)

// Skip unless the world size is at least n (parity-layout tests that scale:
// rank 0 clipper, even ranks mesh 1, odd ranks mesh 2 — valid for any size >= 3
// where both an even and an odd non-zero rank exist).
#define SKIP_UNLESS_MIN_PROCS(clac, n)                                                                                  \
  do {                                                                                                                  \
    const FS_intT nProcs__ = (clac).GetNProcs();                                                                        \
    if(nProcs__ < (n))                                                                                                  \
      GTEST_SKIP() << "needs at least " << (n) << " MPI procs, running with " << nProcs__;                              \
  } while(0)

_FS_END_NAMESPACE
