#include "FSClipping/FSFaceMatcher.h"

#include <gtest/gtest.h>

#include <set>

static FSFaceMatch MakeMatch(const std::vector<FSVec3>& poly)
{
  FSFaceMatch m;
  m.type = FSFaceMatch::INTERSECTING;
  m.clippedPoly3D = poly;
  return m;
}

TEST(AssignGlobalNodeIds, SharedPointsGetSameId)
{
  // Two triangles sharing the edge (1,0,0)-(0,1,0)
  std::vector<FSFaceMatch> matches = {
    MakeMatch({FSVec3(0, 0, 0), FSVec3(1, 0, 0), FSVec3(0, 1, 0)}),
    MakeMatch({FSVec3(1, 0, 0), FSVec3(1, 1, 0), FSVec3(0, 1, 0)}),
  };

  const FS_intT nUnique = FSFaceMatcher::AssignGlobalNodeIds(matches, 1e-12);

  EXPECT_EQ(nUnique, 4);
  ASSERT_EQ(matches[0].nodeGlobalIds.size(), 3u);
  ASSERT_EQ(matches[1].nodeGlobalIds.size(), 3u);

  // (1,0,0): index 1 in match 0, index 0 in match 1
  EXPECT_EQ(matches[0].nodeGlobalIds[1], matches[1].nodeGlobalIds[0]);
  // (0,1,0): index 2 in match 0, index 2 in match 1
  EXPECT_EQ(matches[0].nodeGlobalIds[2], matches[1].nodeGlobalIds[2]);
  // (0,0,0) and (1,1,0) are distinct
  EXPECT_NE(matches[0].nodeGlobalIds[0], matches[1].nodeGlobalIds[1]);
}

TEST(AssignGlobalNodeIds, IdsAreContiguous)
{
  std::vector<FSFaceMatch> matches = {
    MakeMatch({FSVec3(0, 0, 0), FSVec3(1, 0, 0), FSVec3(0, 1, 0)}),
    MakeMatch({FSVec3(1, 0, 0), FSVec3(1, 1, 0), FSVec3(0, 1, 0)}),
    MakeMatch({FSVec3(5, 5, 5), FSVec3(6, 5, 5), FSVec3(5, 6, 5)}),
  };

  const FS_intT nUnique = FSFaceMatcher::AssignGlobalNodeIds(matches, 1e-12);

  // Collect all assigned IDs must be exactly {0, 1, ..., nUnique-1}
  std::set<FS_intT> ids;
  for(const auto& m : matches)
    for(FS_intT id : m.nodeGlobalIds)
      ids.insert(id);

  EXPECT_EQ(static_cast<FS_intT>(ids.size()), nUnique);
  EXPECT_EQ(*ids.begin(), 0);
  EXPECT_EQ(*ids.rbegin(), nUnique - 1);
}

TEST(AssignGlobalNodeIds, ToleranceMergesClosePoints)
{
  const FS_floatT tol = 1e-6;

  // Same point up to 1e-8 (< tol) must be merged.
  std::vector<FSFaceMatch> matches = {
    MakeMatch({FSVec3(0, 0, 0), FSVec3(1, 0, 0), FSVec3(0, 1, 0)}),
    MakeMatch({FSVec3(1e-8, 0, 0), FSVec3(1, 1, 0), FSVec3(0, 1 + 1e-8, 0)}),
  };

  const FS_intT nUnique = FSFaceMatcher::AssignGlobalNodeIds(matches, tol);

  EXPECT_EQ(nUnique, 4);
  EXPECT_EQ(matches[0].nodeGlobalIds[0], matches[1].nodeGlobalIds[0]);
  EXPECT_EQ(matches[0].nodeGlobalIds[2], matches[1].nodeGlobalIds[2]);
}

TEST(AssignGlobalNodeIds, NotCoveredMatchIsNumbered)
{
  FSFaceMatch notCovered;
  notCovered.type = FSFaceMatch::NOT_COVERED;
  notCovered.clippedPoly3D = {FSVec3(0, 0, 0), FSVec3(1, 0, 0), FSVec3(1, 1, 0), FSVec3(0, 1, 0)};

  std::vector<FSFaceMatch> matches = {
    MakeMatch({FSVec3(1, 0, 0), FSVec3(2, 0, 0), FSVec3(1, 1, 0)}),
    notCovered,
  };

  const FS_intT nUnique = FSFaceMatcher::AssignGlobalNodeIds(matches, 1e-12);

  // 5 unique points: (1,0,0) and (1,1,0) shared between the two matches
  EXPECT_EQ(nUnique, 5);
  ASSERT_EQ(matches[1].nodeGlobalIds.size(), 4u);
  EXPECT_EQ(matches[0].nodeGlobalIds[0], matches[1].nodeGlobalIds[1]); // (1,0,0)
  EXPECT_EQ(matches[0].nodeGlobalIds[2], matches[1].nodeGlobalIds[2]); // (1,1,0)
}

TEST(AssignGlobalNodeIds, EmptyMatchesReturnsZero)
{
  std::vector<FSFaceMatch> matches;
  EXPECT_EQ(FSFaceMatcher::AssignGlobalNodeIds(matches, 1e-12), 0);
}
