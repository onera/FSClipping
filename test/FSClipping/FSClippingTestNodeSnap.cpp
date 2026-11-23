#include "FSClac.h"
#include "gtest/gtest.h"

#include "FSClipping/FSClippingUtil.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSNodeSnap.h"

_FS_BEGIN_NAMESPACE

namespace {

// Four unit quads around the centre node c (2x2 patch of [-1,1]^2 in the plane z = 0)
std::vector<FSClippingFace> PatchAround(const FSVec3& c)
{
  const FSVec3 sw(-1, -1, 0), s(0, -1, 0), se(1, -1, 0), w(-1, 0, 0), e(1, 0, 0), nw(-1, 1, 0), n(0, 1, 0),
    ne(1, 1, 0);
  std::vector<FSClippingFace> faces;
  faces.emplace_back(std::vector<FSVec3>{sw, s, c, w});
  faces.emplace_back(std::vector<FSVec3>{s, se, e, c});
  faces.emplace_back(std::vector<FSVec3>{c, e, ne, n});
  faces.emplace_back(std::vector<FSVec3>{w, c, n, nw});
  return faces;
}

FS_floatT SmallestArea(const std::vector<FSFaceMatch>& matches)
{
  FS_floatT smallest = FS_FLOATT_MAX;
  for(const auto& m : matches)
    smallest = std::min(smallest, m.intersectedArea);
  return smallest;
}

} // namespace

const FSVec3 kCentre(0, 0, 0);
const FSVec3 kNearCentre(2e-7, -1.5e-7, 0); // |S R| = 2.5e-7, as on the rotor/stator interface

TEST(FSNodeSnap, NearNodeIsMovedOntoTarget)
{
  const auto moves = FSNodeSnap::Compute(PatchAround(kCentre), PatchAround(kNearCentre), 1e-6);
  ASSERT_EQ(moves.size(), 1u);
  EXPECT_EQ(moves[0].from.Distance(kNearCentre), 0.0);
  EXPECT_EQ(moves[0].to.Distance(kCentre), 0.0);
}

TEST(FSNodeSnap, BelowToleranceOrCoincidentNothingMoves)
{
  EXPECT_TRUE(FSNodeSnap::Compute(PatchAround(kCentre), PatchAround(kNearCentre), 1e-8).empty());
  EXPECT_TRUE(FSNodeSnap::Compute(PatchAround(kCentre), PatchAround(kCentre), 1e-6).empty());
  EXPECT_TRUE(FSNodeSnap::Compute(PatchAround(kCentre), PatchAround(kNearCentre), 0.0).empty());
}

TEST(FSNodeSnap, ApplyMovesFacesAndPoints)
{
  auto moving = PatchAround(kNearCentre);
  const auto moves = FSNodeSnap::Compute(PatchAround(kCentre), moving, 1e-6);
  FSNodeSnap::Apply(moves, moving);
  const auto target = PatchAround(kCentre);
  for(std::size_t f = 0; f < moving.size(); ++f)
    for(std::size_t i = 0; i < moving[f].vertices().size(); ++i)
      EXPECT_EQ(moving[f].vertices()[i].Distance(target[f].vertices()[i]), 0.0);

  std::vector<FSVec3> points{kNearCentre, FSVec3(1, 1, 0)};
  EXPECT_EQ(FSNodeSnap::Apply(moves, points), 1);
  EXPECT_EQ(points[0].Distance(kCentre), 0.0);
  EXPECT_EQ(points[1].Distance(FSVec3(1, 1, 0)), 0.0);
}

// Without snapping the clipping cuts micro-pieces around S and R; after it the two patches are identical
TEST(FSNodeSnap, SnappingRemovesMicroPieces)
{
  FSClac dummyClac;
  const auto subject = PatchAround(kCentre);
  {
    const auto clipped = PatchAround(kNearCentre);
    FSFaceMatcher matcher(dummyClac, subject, clipped, 1e-10);
    std::vector<FSFaceMatch> matches;
    matcher.ComputeMatches(matches);
    EXPECT_GT(matches.size(), 4u);
    EXPECT_LT(SmallestArea(matches), 1e-12);
  }

  auto clipped = PatchAround(kNearCentre);
  FSNodeSnap::Apply(FSNodeSnap::Compute(subject, clipped, 1e-6), clipped);
  FSFaceMatcher matcher(dummyClac, subject, clipped, 1e-10);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);
  ASSERT_EQ(matches.size(), 4u);
  for(const auto& m : matches)
    EXPECT_EQ(m.type, FSFaceMatch::IDENTICAL);
}

_FS_END_NAMESPACE
