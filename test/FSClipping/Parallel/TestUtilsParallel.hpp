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
  mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  return mesh;
}

_FS_END_NAMESPACE
