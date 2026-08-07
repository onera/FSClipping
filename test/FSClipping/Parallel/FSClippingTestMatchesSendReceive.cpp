#include <gtest/gtest.h>
#include "FSClipping/FSClippingInterfacePar.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSClippingFace.h"
#include "TestUtilsParallel.hpp"

// The identity square [-1,1]x[-1,1] is split into four QUAD4 faces.
struct IdentitySquare {
  FSFloatArrayT bottomLeftArray;  // [-1,-1] x [0, 0]
  FSFloatArrayT bottomRightArray; // [ 0,-1] x [1, 0]
  FSFloatArrayT topRightArray;    // [ 0, 0] x [1, 1]
  FSFloatArrayT topLeftArray;     // [-1, 0] x [0, 1]

  IdentitySquare()
  {
    bottomLeftArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT bottomLeft[FSMESH_NNODES_QUAD4][FS_3D] = {{-1, -1, 0}, {0, -1, 0}, {0, 0, 0}, {-1, 0, 0}};

    bottomRightArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT bottomRight[FSMESH_NNODES_QUAD4][FS_3D] = {{0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 0, 0}};

    topRightArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT topRight[FSMESH_NNODES_QUAD4][FS_3D] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};

    topLeftArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT topLeft[FSMESH_NNODES_QUAD4][FS_3D] = {{-1, 0, 0}, {0, 0, 0}, {0, 1, 0}, {-1, 1, 0}};

    for(FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i) {
      for(FS_intT d = 0; d < FS_3D; ++d) {
        bottomLeftArray(i, d) = bottomLeft[i][d];
        bottomRightArray(i, d) = bottomRight[i][d];
        topLeftArray(i, d) = topLeft[i][d];
        topRightArray(i, d) = topRight[i][d];
      }
    }
  }
};

struct ExpectedMatch {
  FS_intT type;
  FS_intT face1;
  FS_intT face2;
};

// Verifies that the received match list corresponds exactly to the expected list.
static void AssertMatches(const std::vector<FSFaceMatch>& matches, const std::vector<ExpectedMatch>& expected)
{
  ASSERT_EQ(matches.size(), expected.size());

  for(std::size_t i = 0; i < matches.size(); ++i) {
    EXPECT_EQ(matches[i].type, expected[i].type);
    EXPECT_EQ(matches[i].face1, expected[i].face1);
  }
}

static void RunParallelMatchTest(FSClac& globalClac, FSClac& localClac, FS_intT meshID,
                                 const std::vector<FSClippingFace>& faces,
                                 const std::vector<ExpectedMatch>& expectedMatches, FS_floatT tol = 1e-12)
{
  // meshID {1,2} send the faces to meshID 0 (procID 0)
  if(meshID == 1 || meshID == 2) {
    FSFaceExchange::GatherSend(globalClac, 0, meshID, faces);
    auto matches = FSMatchExchange::ScatterReceive(globalClac, 0);

    if(meshID == 1)
      AssertMatches(matches, expectedMatches);
    else {
      std::vector<ExpectedMatch> expectedSwapped = expectedMatches;
      for(auto& m : expectedSwapped)
        std::swap(m.face1, m.face2);
      AssertMatches(matches, expectedSwapped);
    }

  }

  // meshID 0 receive the faces and compute the matches
  else {
    auto allGathered = FSFaceExchange::GatherReceiveAll(globalClac);

    FSFaceMatcher matcher(localClac, allGathered.meshA.faces.faces, allGathered.meshB.faces.faces, tol);
    std::vector<FSFaceMatch> matches;
    matcher.ComputeMatches(matches);

    FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshA);
    matcher.ComputeInvertedMatches(matches);
    FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshB);
  }
}

// Resets the face counter before each test
class FSClippingTestMatchesPar : public ::testing::Test
{
protected:
  void SetUp() override { FSClippingFace::ResetIdCounter(); }
};

// Test 1 : Identical
TEST_F(FSClippingTestMatchesPar, Identical)
{
  FSClac globalClac(MPI_COMM_WORLD);

  FS_intT meshID = globalClac.GetProcID();
  FSClac localClac;
  globalClac.DivideIntoGroups(meshID, localClac);

  IdentitySquare sq;
  std::vector<FSClippingFace> faces;

  if(meshID == 1 || meshID == 2) {
    faces.emplace_back(sq.bottomLeftArray);
    faces.emplace_back(sq.bottomRightArray);
    faces.emplace_back(sq.topLeftArray);
    faces.emplace_back(sq.topRightArray);
  }

  // Each subject face i is identical to clipped face i
  const std::vector<ExpectedMatch> expected = {
    {FSFaceMatch::IDENTICAL, 0, 0},
    {FSFaceMatch::IDENTICAL, 1, 1},
    {FSFaceMatch::IDENTICAL, 2, 2},
    {FSFaceMatch::IDENTICAL, 3, 3},
  };

  RunParallelMatchTest(globalClac, localClac, meshID, faces, expected);
}

// Test 2: intersection
TEST_F(FSClippingTestMatchesPar, Intersection)
{
  FSClac globalClac(MPI_COMM_WORLD);

  FS_intT meshID = globalClac.GetProcID();
  FSClac localClac;
  globalClac.DivideIntoGroups(meshID, localClac);

  IdentitySquare sq;
  std::vector<FSClippingFace> faces;

  FSFloatArrayT fullArray(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT fullSquare[FSMESH_NNODES_QUAD4][FS_3D] = {{-1.0, -1.0, 0}, {1.0, -1.0, 0}, {1.0, 1.0, 0}, {-1.0, 1.0, 0}};

  if(meshID == 1) {
    faces.emplace_back(sq.bottomLeftArray);  // id 0
    faces.emplace_back(sq.bottomRightArray); // id 1
    faces.emplace_back(sq.topLeftArray);     // id 2
    faces.emplace_back(sq.topRightArray);    // id 3
  } else if(meshID == 2) {
    for(FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i)
      for(FS_intT d = 0; d < FS_3D; ++d)
        fullArray(i, d) = fullSquare[i][d];

    faces.emplace_back(fullArray); // id 0
  }

  // Each subject quadrant intersects the single clipped face
  const std::vector<ExpectedMatch> expected = {
    {FSFaceMatch::INTERSECTING, 0, 0}, // bottomLeft
    {FSFaceMatch::INTERSECTING, 1, 0}, // bottomRight
    {FSFaceMatch::INTERSECTING, 2, 0}, // topLeft
    {FSFaceMatch::INTERSECTING, 3, 0}, // topRight
  };

  RunParallelMatchTest(globalClac, localClac, meshID, faces, expected);
}