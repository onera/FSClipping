// FSClippingTestMatchesPar.cpp

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
    FS_floatT bottomLeft[FSMESH_NNODES_QUAD4][FS_3D] = {
      {-1, -1, 0}, {0, -1, 0}, {0, 0, 0}, {-1, 0, 0}};

    bottomRightArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT bottomRight[FSMESH_NNODES_QUAD4][FS_3D] = {
      {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 0, 0}};

    topRightArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT topRight[FSMESH_NNODES_QUAD4][FS_3D] = {
      {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};

    topLeftArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT topLeft[FSMESH_NNODES_QUAD4][FS_3D] = {
      {-1, 0, 0}, {0, 0, 0}, {0, 1, 0}, {-1, 1, 0}};

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


// ─── Expected match helper ────────────────────────────────────────────────────
// Describes the full expected state of one FSFaceMatch.
// face1 = index of the subject face in the subject vector
// face2 = index of the clipped face in the clipped vector
struct ExpectedMatch {
  FS_intT type;
  FS_intT face1;
  FS_intT face2;
};

// Verifies that the received match list corresponds exactly to the expected list.
// Matches are sorted by (face1, face2) before comparison so the check
// is order-independent.
static void AssertMatches(const std::vector<FSFaceMatch>& matches,
                          const std::vector<ExpectedMatch>& expected)
{
  ASSERT_EQ(matches.size(), expected.size());

  // Sort received matches by (owner1, owner2) for deterministic comparison
  // std::vector<FSFaceMatch> sorted = matches;
  // std::sort(sorted.begin(), sorted.end(),
  //           [](const FSFaceMatch& a, const FSFaceMatch& b) {
  //             return a.face1 != b.face1
  //                      ? a.face1 < b.face1
  //                      : a.face2 < b.face2;
  //           });

  // // Sort expected the same way
  // std::vector<ExpectedMatch> sortedExp = expected;
  // std::sort(sortedExp.begin(), sortedExp.end(),
  //           [](const ExpectedMatch& a, const ExpectedMatch& b) {
  //             return a.face1 != b.face1
  //                      ? a.face1 < b.face1
  //                      : a.face2 < b.face2;
  //           });

  for(std::size_t i = 0; i < matches.size(); ++i) {
    EXPECT_EQ(matches[i].type, expected[i].type);

    EXPECT_EQ(matches[i].face1, expected[i].face1);

    // EXPECT_EQ(sorted[i].face2, sortedExp[i].face2)
    //  << "Match " << i << ": wrong face2";
  }
}

// ─── Updated parallel helper ──────────────────────────────────────────────────
// Same 3-process protocol as before, but now takes the full expected match list
// instead of just a type + count.
static void RunParallelMatchTest(FSClac& globalClac,
                                 FSClac& localClac,
                                 const std::vector<FSClippingFace>& subjectFaces,
                                 const std::vector<FSClippingFace>& clippedFaces,
                                 const std::vector<ExpectedMatch>& expectedMatches,
                                 FS_floatT tol = 1e-12)
{
  const FS_intT rank = globalClac.GetProcID();

  // ── Proc 0: send subject faces, receive and assert matches ───────────────
  if(rank == 0) {
    FSFaceExchange::Send(globalClac, 2, subjectFaces);
    std::vector<FSFaceMatch> matches = FSMatchExchange::Receive(globalClac, 2);
    AssertMatches(matches, expectedMatches);
  }

  // ── Proc 1: send clipped faces, receive and assert matches ───────────────
  else if(rank == 1) {
    FSFaceExchange::Send(globalClac, 2, clippedFaces);
    std::vector<FSFaceMatch> matches = FSMatchExchange::Receive(globalClac, 2);
    AssertMatches(matches, expectedMatches);
  }

  // ── Proc 2: compute matches and broadcast to Proc 0 and Proc 1 ──────────
  else {
    auto subjectReceived = FSFaceExchange::Receive(globalClac, 0);
    auto clippedReceived = FSFaceExchange::Receive(globalClac, 1);

    FSFaceMatcher matcher(localClac,
                          subjectReceived.faces,
                          clippedReceived.faces,
                          tol);
    std::vector<FSFaceMatch> matches;
    matcher.ComputeMatches(matches);

    FSMatchExchange::Send(globalClac, 0, matches);
    FSMatchExchange::Send(globalClac, 1, matches);
  }
}

// FSClippingFace's coordinate-only constructors number faces with a process-wide
// static counter (s_nextId_) that keeps incrementing across TEST bodies run in
// the same process. These tests assert on absolute faceIndex values, so the
// fixture resets the counter before each test — making them independent of run
// order and of whether they run one-per-process (CTest) or all in one process.
class FSClippingTestMatchesPar : public ::testing::Test
{
protected:
  void SetUp() override { FSClippingFace::ResetIdCounter(); }
};

// ─── Test 1 : Identical(parallel) ────────────────────────────────────────────
// Subject and clipped are the same four faces inserted in the same order.
// Each face i matches itself → face1 == face2 == i for all 4 matches.
//
// Subject vector layout:
//   index 0 → bottomLeft   [-1,-1]x[0,0]
//   index 1 → bottomRight  [ 0,-1]x[1,0]
//   index 2 → topLeft      [-1, 0]x[0,1]
//   index 3 → topRight     [ 0, 0]x[1,1]
TEST_F(FSClippingTestMatchesPar, Identical)
{
  FSClac globalClac(MPI_COMM_WORLD);

  FS_intT meshID = globalClac.GetProcID();
  FSClac localClac;
  globalClac.DivideIntoGroups(meshID, localClac);

  IdentitySquare sq;

  std::vector<FSClippingFace> subjectFaces;
  subjectFaces.emplace_back(sq.bottomLeftArray);  // owner 0
  subjectFaces.emplace_back(sq.bottomRightArray); // owner 1
  subjectFaces.emplace_back(sq.topLeftArray);     // owner 2
  subjectFaces.emplace_back(sq.topRightArray);    // owner 3

  std::vector<FSClippingFace> clippedFaces = subjectFaces;

  // Each subject face i is identical to clipped face i
  const std::vector<ExpectedMatch> expected = {
    {FSFaceMatch::IDENTICAL, 1, 1},
    {FSFaceMatch::IDENTICAL, 2, 2},
    {FSFaceMatch::IDENTICAL, 3, 3},
    {FSFaceMatch::IDENTICAL, 4, 4},
  };

  RunParallelMatchTest(globalClac, localClac,
                       subjectFaces, clippedFaces, expected);
}

// ─── Test 2: Intersection (parallel) ─────────────────────────────────────────
// One clipped quad (index 0) crosses all four quadrants of the subject mesh.
// Each subject face (index 0..3) intersects the single clipped face → owner2 is
// always 0 for every match.
//
// Subject vector layout: same as Test 1 (index 0..3).
// Clipped vector layout:
//   index 0 → centerSquare  [-0.5,-0.5]x[0.5,1]
TEST_F(FSClippingTestMatchesPar, Intersection)
{

  FSClac globalClac(MPI_COMM_WORLD);

  FS_intT meshID = globalClac.GetProcID();
  FSClac localClac;
  globalClac.DivideIntoGroups(meshID, localClac);

  IdentitySquare sq;

  std::vector<FSClippingFace> subjectFaces;
  subjectFaces.emplace_back(sq.bottomLeftArray);  // owner 5
  subjectFaces.emplace_back(sq.bottomRightArray); // owner 6
  subjectFaces.emplace_back(sq.topLeftArray);     // owner 7
  subjectFaces.emplace_back(sq.topRightArray);    // owner 8

  FSFloatArrayT centerArray(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT centerSquare[FSMESH_NNODES_QUAD4][FS_3D] = {
    {-0.5, -0.5, 0}, {0.5, -0.5, 0}, {0.5, 0.5, 0}, {-0.5, 1, 0}};
  for(FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i)
    for(FS_intT d = 0; d < FS_3D; ++d)
      centerArray(i, d) = centerSquare[i][d];

  std::vector<FSClippingFace> clippedFaces;
  clippedFaces.emplace_back(centerArray); // owner 9

  // The single clipped face intersects every subject quadrant, but covers
  // none of them entirely: each subject face is kept unclipped (NOT_COVERED).
  const std::vector<ExpectedMatch> expected = {
    {FSFaceMatch::NOT_COVERED, 1, -1}, // bottomLeft  partially covered
    {FSFaceMatch::NOT_COVERED, 2, -1}, // bottomRight partially covered
    {FSFaceMatch::NOT_COVERED, 3, -1}, // topLeft     partially covered
    {FSFaceMatch::NOT_COVERED, 4, -1}, // topRight    partially covered
  };

  RunParallelMatchTest(globalClac, localClac,
                       subjectFaces, clippedFaces, expected);
}

// ─── Test 3: Included (parallel) ─────────────────────────────────────────────
// One small clipped quad (index 0) lies entirely within subject face 0
// (bottomLeft). No other subject face overlaps it, so exactly one match is
// produced: face1 == 0, face2 == 0.
//
// Subject vector layout: same as Test 1 (index 0..3).
// Clipped vector layout:
//   index 0 → includedSquare  [-0.75,-0.75]x[-0.25,-0.25]
TEST_F(FSClippingTestMatchesPar, Included)
{
  FSClac globalClac(MPI_COMM_WORLD);

  FS_intT meshID = globalClac.GetProcID();
  FSClac localClac;
  globalClac.DivideIntoGroups(meshID, localClac);


  IdentitySquare sq;

  std::vector<FSClippingFace> subjectFaces;
  subjectFaces.emplace_back(sq.bottomLeftArray);  // owner 1
  subjectFaces.emplace_back(sq.bottomRightArray); // owner 2
  subjectFaces.emplace_back(sq.topLeftArray);     // owner 3
  subjectFaces.emplace_back(sq.topRightArray);    // owner 4

  FSFloatArrayT includedArray(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT includedSquare[FSMESH_NNODES_QUAD4][FS_3D] = {
    {-0.75, -0.75, 0}, {-0.25, -0.75, 0}, {-0.25, -0.25, 0}, {-0.75, -0.25, 0}};
  for(FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i)
    for(FS_intT d = 0; d < FS_3D; ++d)
      includedArray(i, d) = includedSquare[i][d];

  std::vector<FSClippingFace> clippedFaces;
  clippedFaces.emplace_back(includedArray); // owner 14

  // No subject face is fully covered: bottomLeft only partially (the small
  // quad lies inside it) and the three others not at all. Every subject face
  // is therefore kept unclipped as a NOT_COVERED match.
  const std::vector<ExpectedMatch> expected = {
    {FSFaceMatch::NOT_COVERED, 1, -1},
    {FSFaceMatch::NOT_COVERED, 2, -1},
    {FSFaceMatch::NOT_COVERED, 3, -1},
    {FSFaceMatch::NOT_COVERED, 4, -1},
  };

  RunParallelMatchTest(globalClac, localClac,
                       subjectFaces, clippedFaces, expected, 1e-6);
}
