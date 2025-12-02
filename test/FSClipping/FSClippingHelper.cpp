#include "FSBVHTree.h"
#include "FSClac.h"
#include "FSClipping/FSClippingFace.h"
#include "FSCommon.h"
#include "FSError.h"
#include "FSLog.h"

#include "FSMesh/FSBoundingBoxUtil.h"
#include "FSMesh/FSHashableFace.h"
#include "FSMesh/FSMesh.h"
#include "FSMesh/FSMeshCheck.h"
#include "FSMesh/FSMeshData.h"
#include "FSMesh/FSMeshFaceExtractor.h"
#include "FSMesh/FSMeshImportParamsTAU.h"
#include "FSMesh/FSMeshOpParams.h"
#include "FSMesh/FSMeshOpParamsArray.h"
#include "FSMesh/FSMeshPrintInfo.h"

#include "gtest/gtest.h"
#include <ostream>

#include "FSGeometry/FSPlaneSurface.h"

#include "FSClipping/FSBoundaryFace.h"
#include "FSClipping/FSCellAdress.h"
#include "FSClipping/FSClippingAlgo.h"
#include "FSClipping/FSClippingUtil.h"
#include "FSClipping/FSFace.h"
#include "FSClipping/FSFaceSeparator.h"

_FS_BEGIN_NAMESPACE

TEST(FSBoundingBoxIntersec, IntersectingBox) {
  FSClac clac;
  FSMeshImportParamsTAU params;
  params.mMeshFilename = "${HOME}/test/FSDM_script/hexa.grid";

  FSMesh mesh(&clac);
  ASSERT_TRUE(mesh.ImportMesh(&params));

  FSMeshFaceExtractor ex;
  mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  ASSERT_TRUE(ex.PrepareFaceConnectivity(mesh.GetMeshData()->GetUnstructCells(),
                                         true, false));
  ASSERT_TRUE(ex.PrepareFaceNodeCoordinates(FSQuantityDescArrayT()));

  auto bdryFaces =
      FSFaceSeparator::SeparateBoundariesFaceWithMarker(mesh, ex, 1);
  ASSERT_FALSE(bdryFaces.empty());

  std::vector<FSClippingFace> clipFaces;

  for (const auto &f : bdryFaces) {
    clipFaces.emplace_back(f, ex);

    // tests géométriques
    EXPECT_NEAR(clipFaces.back().tangentVector1().L2Norm(), 1.0, 1e-12);
    EXPECT_NEAR(clipFaces.back().normal().L2Norm(), 1.0, 1e-12);

    //  const auto &fBoundingBox = clipFaces.back().boundingBox();
    //  std::cout << ex.GetFaceIndex(*(f._faceFSDM)) << " : "
    //            << fBoundingBox.boxMinMax[0] << ", " <<
    //            fBoundingBox.boxMinMax[1]
    //            << ", " << fBoundingBox.boxMinMax[2] << std::endl;
  }

  FSFloatArrayT boundingBoxes(clipFaces.size(), 6);
  FSIntArrayT localKeys(clipFaces.size());
  for (std::size_t i = 0; i < clipFaces.size(); i++) {
    boundingBoxes(i, 0) = clipFaces[i].boundingBox().boxMinMax[0];
    boundingBoxes(i, 1) = clipFaces[i].boundingBox().boxMinMax[1];
    boundingBoxes(i, 2) = clipFaces[i].boundingBox().boxMinMax[2];
    boundingBoxes(i, 3) = clipFaces[i].boundingBox().boxMinMax[3];
    boundingBoxes(i, 4) = clipFaces[i].boundingBox().boxMinMax[4];
    boundingBoxes(i, 5) = clipFaces[i].boundingBox().boxMinMax[5];
    localKeys(i) = ex.GetFaceIndex(*(clipFaces[i].topo()._faceFSDM));
  }

  FSBVHTree facebvhtree;
  FSBoundingBoxUtil::GatherBoundingBoxesIntoBVHTree(&clac, boundingBoxes,
                                                    localKeys, facebvhtree);

  FS_floatT box[2 * FS_3D] = {-1, -0.75, -0.75, -1, -0.25, -0.25};
  FSIntArrayT outIndices;

  const FS_intT numFound =
      facebvhtree.FindBoxesIntersectingWithBox(box, outIndices);

  FSIntArrayT outKeys = facebvhtree.GetKeysForCells(outIndices);

  ASSERT_TRUE(outKeys.Size() == 4);

  FSFloatArrayT coord;
  // for (FS_intT i = 0; i < outKeys.Size(); i++) {
  //   //ex.GetFaceNodeCoordinates(outKeys[i], coord);
  //   std::cout << outKeys[i] << std::endl;
  // }
}
_FS_END_NAMESPACE
