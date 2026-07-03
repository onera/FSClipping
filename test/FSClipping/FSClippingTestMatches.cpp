#include "FSClac.h"
#include "FSCommon.h"
#include "FSError.h"
#include "FSLog.h"
#include "FSPlaneSurface.h"

#include "gtest/gtest.h"

#include "FSClipping/FSClippingFace.h"
#include "FSClipping/FSFaceMatcher.h"

_FS_BEGIN_NAMESPACE

struct IdentitySquare {
  // square [-1,-1]x[0,0]
  FSFloatArrayT bottomLeftArray;
  // square [0,-1]x[1,0]
  FSFloatArrayT bottomRightArray;
  // square [0,0]x[1,1]
  FSFloatArrayT topRightArray;
  // square [-1,0]x[0,1]
  FSFloatArrayT topLeftArray;

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

struct IdentitySquareZ {
  // square [-1,-1]x[0,0]
  FSFloatArrayT bottomLeftArray;
  // square [0,-1]x[1,0]
  FSFloatArrayT bottomRightArray;
  // square [0,0]x[1,1]
  FSFloatArrayT topRightArray;
  // square [-1,0]x[0,1]
  FSFloatArrayT topLeftArray;

  IdentitySquareZ()
  {

    bottomLeftArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT bottomLeft[FSMESH_NNODES_QUAD4][FS_3D] = {{0, -1, -1}, {0, -1, 0}, {0, 0, 0}, {0, 0, -1}};

    bottomRightArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT bottomRight[FSMESH_NNODES_QUAD4][FS_3D] = {{0, -1, 0}, {0, -1, 1}, {0, 0, 1}, {0, 0, 0}};

    topRightArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT topRight[FSMESH_NNODES_QUAD4][FS_3D] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}};

    topLeftArray.Resize(FSMESH_NNODES_QUAD4, FS_3D);
    FS_floatT topLeft[FSMESH_NNODES_QUAD4][FS_3D] = {{0, 0, -1}, {0, 0, 0}, {0, 1, 0}, {0, 1, -1}};

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

// --- Test 1: identical
TEST(FSClippingTestMatches, Identical)
{

  FSClac dummyClac;
  IdentitySquare identitySquare;
  std::vector<FSClippingFace> subject;

  subject.emplace_back(identitySquare.bottomLeftArray);
  subject.emplace_back(identitySquare.bottomRightArray);
  subject.emplace_back(identitySquare.topLeftArray);
  subject.emplace_back(identitySquare.topRightArray);

  std::vector<FSClippingFace> clipped = subject;

  // Build matcher
  FSFaceMatcher matcherIdentical(dummyClac, subject, clipped);
  std::vector<FSFaceMatch> matches;
  matcherIdentical.ComputeMatches(matches);
  ASSERT_TRUE(matches.size() == 4);
  for(const auto& m : matches)
    ASSERT_TRUE(m.type == FSFaceMatch::IDENTICAL);
}

// --- Test 2: intersection
TEST(FSClippingTestMatches, Intersection)
{

  FSClac dummyClac;
  IdentitySquare identitySquare;
  std::vector<FSClippingFace> subject;
  std::vector<FSClippingFace> clipped;

  FSFloatArrayT centerSquareArray(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT centerSquare[FSMESH_NNODES_QUAD4][FS_3D] = {{-0.5, -0.5, 0}, {0.5, -0.5, 0}, {0.5, 0.5, 0}, {-0.5, 1, 0}};

  for(FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i)
    for(FS_intT d = 0; d < FS_3D; ++d)
      centerSquareArray(i, d) = centerSquare[i][d];

  clipped.emplace_back(centerSquareArray);
  subject.emplace_back(identitySquare.bottomLeftArray);
  subject.emplace_back(identitySquare.bottomRightArray);
  subject.emplace_back(identitySquare.topLeftArray);
  subject.emplace_back(identitySquare.topRightArray);

  // Build matcher
  FSFaceMatcher matcherIntersection(dummyClac, subject, clipped);
  std::vector<FSFaceMatch> matches;
  matcherIntersection.ComputeMatches(matches);
  ASSERT_TRUE(matches.size() == 4);

  // Each subject quadrant is only partially covered by the clipped face, so
  // the clipping is discarded and the original face is kept (NOT_COVERED).
  for(const auto& m : matches) {
    ASSERT_TRUE(m.type == FSFaceMatch::NOT_COVERED);
    ASSERT_TRUE(m.face2 == -1);
    ASSERT_TRUE(m.clippedPoly3D.size() == 4);
  }
}

// --- Test 2: intersection along z
TEST(FSClippingTestMatches, IntersectionZ)
{

  FSClac dummyClac;
  IdentitySquareZ identitySquare;
  std::vector<FSClippingFace> subject;
  std::vector<FSClippingFace> clipped;

  FSFloatArrayT centerSquareArray(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT centerSquare[FSMESH_NNODES_QUAD4][FS_3D] = {{0, -0.5, -0.5}, {0, -0.5, 0.5}, {0, 0.5, 0.5}, {0, 1, -0.5}};

  for(FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i)
    for(FS_intT d = 0; d < FS_3D; ++d)
      centerSquareArray(i, d) = centerSquare[i][d];

  clipped.emplace_back(centerSquareArray);
  subject.emplace_back(identitySquare.bottomLeftArray);
  subject.emplace_back(identitySquare.bottomRightArray);
  subject.emplace_back(identitySquare.topLeftArray);
  subject.emplace_back(identitySquare.topRightArray);

  // Build matcher
  FSFaceMatcher matcherIntersection(dummyClac, subject, clipped);
  std::vector<FSFaceMatch> matches;
  matcherIntersection.ComputeMatches(matches);
  ASSERT_TRUE(matches.size() == 4);

  // Each subject quadrant is only partially covered by the clipped face, so
  // the clipping is discarded and the original face is kept (NOT_COVERED).
  for(const auto& m : matches) {
    ASSERT_TRUE(m.type == FSFaceMatch::NOT_COVERED);
    ASSERT_TRUE(m.face2 == -1);
    ASSERT_TRUE(m.clippedPoly3D.size() == 4);
  }
}

// --- Test: inverted matches (coverage is decided independently on each side)
TEST(FSClippingTestMatches, InvertedMatches)
{

  FSClac dummyClac;
  IdentitySquare identitySquare;
  std::vector<FSClippingFace> subject;
  std::vector<FSClippingFace> clipped;

  // Center square [-0.5,0.5]x[-0.5,0.5]: crosses the four subject quadrants
  // and is entirely tiled by its four intersections with them.
  FSFloatArrayT centerSquareArray(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT centerSquare[FSMESH_NNODES_QUAD4][FS_3D] = {
    {-0.5, -0.5, 0}, {0.5, -0.5, 0}, {0.5, 0.5, 0}, {-0.5, 0.5, 0}};

  for(FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i)
    for(FS_intT d = 0; d < FS_3D; ++d)
      centerSquareArray(i, d) = centerSquare[i][d];

  clipped.emplace_back(centerSquareArray);
  subject.emplace_back(identitySquare.bottomLeftArray);
  subject.emplace_back(identitySquare.bottomRightArray);
  subject.emplace_back(identitySquare.topLeftArray);
  subject.emplace_back(identitySquare.topRightArray);

  FSFaceMatcher matcher(dummyClac, subject, clipped);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  // Subject side: each quadrant is only partially covered by the center
  // square, so the original faces are kept (NOT_COVERED).
  ASSERT_TRUE(matches.size() == 4);
  for(const auto& m : matches)
    ASSERT_TRUE(m.type == FSFaceMatch::NOT_COVERED);

  // Clipped side: the center square is fully covered by its four
  // intersections, so the raw INTERSECTING matches are kept.
  std::vector<FSFaceMatch> inverted;
  matcher.ComputeInvertedMatches(inverted);
  ASSERT_TRUE(inverted.size() == 4);
  for(const auto& m : inverted) {
    ASSERT_TRUE(m.type == FSFaceMatch::INTERSECTING);
    ASSERT_TRUE(m.face1 == clipped[0].faceIndex());
    ASSERT_TRUE(m.clippedPoly3D.size() >= 3);
  }
}

// --- Test 2: included
TEST(FSClippingTestMatches, included)
{

  FSClac dummyClac;
  IdentitySquare identitySquare;
  std::vector<FSClippingFace> subject;
  std::vector<FSClippingFace> clipped;

  FSFloatArrayT includedSquareArray(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT includedSquare[FSMESH_NNODES_QUAD4][FS_3D] = {
    {-0.75, -0.75, 0}, {-0.25, -0.75, 0}, {-0.25, -0.25, 0}, {-0.75, -0.25, 0}};
  // should be included in the bottom left square of the identity square

  for(FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i)
    for(FS_intT d = 0; d < FS_3D; ++d)
      includedSquareArray(i, d) = includedSquare[i][d];

  clipped.emplace_back(includedSquareArray);
  subject.emplace_back(identitySquare.bottomLeftArray);
  subject.emplace_back(identitySquare.bottomRightArray);
  subject.emplace_back(identitySquare.topLeftArray);
  subject.emplace_back(identitySquare.topRightArray);

  // Build matcher
  FSFaceMatcher matcherIntersection(dummyClac, subject, clipped, 1e-6);
  std::vector<FSFaceMatch> matches;
  matcherIntersection.ComputeMatches(matches);

  // No subject face is fully covered: bottomLeft only partially (the small
  // quad lies inside it) and the three others not at all. Every subject face
  // is therefore kept unclipped as a NOT_COVERED match.
  ASSERT_TRUE(matches.size() == 4);

  for(const auto& m : matches) {
    ASSERT_TRUE(m.type == FSFaceMatch::NOT_COVERED);
    ASSERT_TRUE(m.face2 == -1);
    ASSERT_TRUE(m.clippedPoly3D.size() == 4);
  }
}

_FS_END_NAMESPACE
