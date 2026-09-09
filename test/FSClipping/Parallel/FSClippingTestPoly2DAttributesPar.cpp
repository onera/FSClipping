// Attributes of the clipped Poly2D cells, run at exactly 9 MPI procs:
// rank 0 = clipper, ranks 1..8 split between the two meshes. 9 procs is the smallest count at which RCB splits a
// volume/surface pair on these meshes, which is what this test is about.
//
// A Poly2D is born from a match, and its parent is the surface cell the matched face closes. When RCB splits a
// volume/surface pair, that parent is owned by another proc: its index means nothing in the local pool, so the
// CADGroupID cannot be read back locally and has to travel with the match. Without that, the copy either reads an
// unrelated local cell or runs past the end of the pool.

#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"
#include <FSDataManagerData.h>

_FS_BEGIN_NAMESPACE

namespace {

/// Counts the Poly2D cells of `mesh` whose CADGroupID differs from `expectedMarker`.
FS_int32T CountWrongPoly2DMarkers(const FSMesh& mesh, FS_intT expectedMarker)
{
  static const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  const FSUnstructMeshData& meshData = mesh.GetMeshData()->GetUnstructCells();

  const FSCellPool* cellPool = meshData.GetCellPool(FSMeshEnums::CT_Poly2D);
  if(!cellPool || meshData.GetCellAttributes(FSMeshEnums::CT_Poly2D).Size() == 0)
    return 0;

  const FSIntArrayT& markers = meshData.GetCellAttribute(attrCADGroupIDName, FSMeshEnums::CT_Poly2D);

  FS_int32T wrong = 0;
  for(FS_intT i = 0; i < markers.Size(); ++i)
    if(markers[i] != expectedMarker)
      ++wrong;
  return wrong;
}

} // namespace

TEST(FSClippingTestPoly2DAttributesPar, ClippedFacesKeepTheirCADGroupID)
{
  FSClac globalClac(MPI_COMM_WORLD);
  SKIP_UNLESS_EXACT_PROCS(globalClac, 9);

  const FS_intT procId = globalClac.GetProcID();
  const FS_intT meshID = SymMeshID(procId);

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSDataManagerData* data = new FSDataManagerData(&globalClac);

  FSClippedMeshParams clippedParams;
  clippedParams.mMeshKeyOriginal1 = "original1";
  clippedParams.mMeshKeyOriginal2 = "original2";
  clippedParams.mMeshKeyClipped1 = "clippedMesh1";
  clippedParams.mMeshKeyClipped2 = "clippedMesh2";
  clippedParams.mMarker1 = 6;
  clippedParams.mMarker2 = 5;
  clippedParams.mTol = 1e-8;

  if(meshID == 1) {
    FSMesh* meshPtr = data->GetMesh("original1", &clac);
    LoadMeshWithClac(meshPtr, MeshPath("input/cube_hexa_coarse_par.grid"));
    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  }

  if(meshID == 2) {
    FSMesh* meshPtr = data->GetMesh("original2", &clac);
    LoadMeshWithClac(meshPtr, MeshPath("input/cube_hexa_fine_par.grid"));
    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  }

  FSClippedMesh clippedMesh(&globalClac);
  ASSERT_TRUE(clippedMesh.DoOp(data, &clippedParams));

  if(meshID == 0)
    return; // clipper proc holds no clipped mesh

  const FSString clippedKey = (meshID == 1) ? "clippedMesh1" : "clippedMesh2";
  const FS_intT expectedMarker = (meshID == 1) ? clippedParams.mMarker1 : clippedParams.mMarker2;

  FS_int32T localWrong = CountWrongPoly2DMarkers(*data->GetMesh(clippedKey, &clac), expectedMarker);
  FS_int32T totalWrong = 0;
  clac.Sum(&localWrong, &totalWrong, 1);

  // Every clipped face lies on the clipped boundary, so all Poly2D must carry that marker -- those whose parent
  // surface cell was left on another proc included.
  EXPECT_EQ(totalWrong, 0);
}

_FS_END_NAMESPACE
