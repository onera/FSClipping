#include "FSClac.h"
#include "gtest/gtest.h"

#include "FSClipping/FSFaceMatcher.h"

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

FS_floatT AreaXY(const std::vector<FSVec3>& poly)
{
  FS_floatT a = 0.0;
  for(std::size_t i = 0; i < poly.size(); ++i) {
    const FSVec3& p = poly[i];
    const FSVec3& q = poly[(i + 1) % poly.size()];
    a += p[0] * q[1] - q[0] * p[1];
  }
  return 0.5 * std::abs(a);
}

bool IsVertexOf(const FSVec3& p, const std::vector<FSClippingFace>& faces)
{
  for(const auto& f : faces)
    for(const auto& v : f.vertices())
      if(v.Distance(p) == 0.0)
        return true;
  return false;
}

} // namespace

const FSVec3 kCentre(0, 0, 0);
const FSVec3 kNearCentre(2e-7, -1.5e-7, 0); // |S R| = 2.5e-7, as on the rotor/stator interface

TEST(FSVertexSnap, DisabledKeepsMicroPieces)
{
  FSClac dummyClac;
  const auto subject = PatchAround(kCentre), clipped = PatchAround(kNearCentre);
  FSFaceMatcher matcher(dummyClac, subject, clipped, 1e-10);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);
  EXPECT_GT(matches.size(), 4u);
  EXPECT_EQ(matcher.SnapDropped(0), 0);
}

// Node-node: each side keeps its own node, S for the subject, R for the clipped mesh
TEST(FSVertexSnap, NearNodesMergeOnTheirOwnNode)
{
  FSClac dummyClac;
  const auto subject = PatchAround(kCentre), clipped = PatchAround(kNearCentre);
  FSFaceMatcher matcher(dummyClac, subject, clipped, 1e-10, 1e-6);

  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);
  ASSERT_EQ(matches.size(), 4u);
  for(const auto& m : matches) {
    EXPECT_GT(AreaXY(m.clippedPoly3D), 0.9);
    for(const auto& p : m.clippedPoly3D)
      EXPECT_TRUE(IsVertexOf(p, subject));
  }

  std::vector<FSFaceMatch> inverted;
  matcher.ComputeInvertedMatches(inverted);
  ASSERT_EQ(inverted.size(), 4u);
  for(const auto& m : inverted) {
    EXPECT_GT(AreaXY(m.clippedPoly3D), 0.9);
    for(const auto& p : m.clippedPoly3D)
      EXPECT_TRUE(IsVertexOf(p, clipped));
  }
  EXPECT_GT(matcher.SnapDropped(0), 0);
  EXPECT_EQ(matcher.SnapDropped(0), matcher.SnapDropped(1));
}

// Node-edge: R lies 2e-7 off the subject edge x = 0, no subject node is near it
TEST(FSVertexSnap, NodeNearEdgeMergesOnTheEdge)
{
  const FSVec3 R(2e-7, 0.3, 0);
  std::vector<FSClippingFace> subject, clipped;
  subject.emplace_back(std::vector<FSVec3>{{-1, -1, 0}, {0, -1, 0}, {0, 1, 0}, {-1, 1, 0}});
  subject.emplace_back(std::vector<FSVec3>{{0, -1, 0}, {1, -1, 0}, {1, 1, 0}, {0, 1, 0}});
  clipped.emplace_back(std::vector<FSVec3>{{-1, -1, 0}, {0, -1, 0}, R, {-1, 0.3, 0}});
  clipped.emplace_back(std::vector<FSVec3>{{0, -1, 0}, {1, -1, 0}, {1, 0.3, 0}, R});
  clipped.emplace_back(std::vector<FSVec3>{R, {1, 0.3, 0}, {1, 1, 0}, {0, 1, 0}});
  clipped.emplace_back(std::vector<FSVec3>{{-1, 0.3, 0}, R, {0, 1, 0}, {-1, 1, 0}});

  FSClac dummyClac;
  FSFaceMatcher matcher(dummyClac, subject, clipped, 1e-10, 1e-6);
  std::vector<FSFaceMatch> matches, inverted;
  matcher.ComputeMatches(matches);
  matcher.ComputeInvertedMatches(inverted);

  ASSERT_EQ(matches.size(), 4u);
  ASSERT_EQ(inverted.size(), 4u);
  bool clippedSeesR = false;
  for(const auto& m : matches) {
    EXPECT_GT(AreaXY(m.clippedPoly3D), 0.25);
    for(const auto& p : m.clippedPoly3D) // the subject side stays on its own edge
      EXPECT_TRUE(std::abs(p[0]) < 1e-14 || std::abs(p[0]) > 1e-3) << p[0];
  }
  for(const auto& m : inverted) {
    EXPECT_GT(AreaXY(m.clippedPoly3D), 0.25);
    for(const auto& p : m.clippedPoly3D)
      clippedSeesR |= (p.Distance(R) == 0.0);
  }
  EXPECT_TRUE(clippedSeesR);
}

// Refined zone (rotor/stator case, edges of 2e-6): each rotor node of a thin strip is 2e-7 off its stator node,
// the two rotor nodes are 1.6e-6 apart, below snapTol: each pairs with its own stator node, never with the other
TEST(FSVertexSnap, ThinStripKeepsOneNodePerGroup)
{
  const FS_floatT w = 2e-6, off = 2e-7;
  auto strips = [](FS_floatT a, FS_floatT b) {
    std::vector<FSClippingFace> faces;
    faces.emplace_back(std::vector<FSVec3>{{-1, -1, 0}, {a, -1, 0}, {a, 1, 0}, {-1, 1, 0}});
    faces.emplace_back(std::vector<FSVec3>{{a, -1, 0}, {b, -1, 0}, {b, 1, 0}, {a, 1, 0}});
    faces.emplace_back(std::vector<FSVec3>{{b, -1, 0}, {1, -1, 0}, {1, 1, 0}, {b, 1, 0}});
    return faces;
  };
  const auto subject = strips(0.0, w), clipped = strips(off, w - off);

  FSClac dummyClac;
  FSFaceMatcher matcher(dummyClac, subject, clipped, 1e-10, 1.7e-6);
  std::vector<FSFaceMatch> matches, inverted;
  matcher.ComputeMatches(matches);
  matcher.ComputeInvertedMatches(inverted);

  ASSERT_EQ(matches.size(), 3u);
  ASSERT_EQ(inverted.size(), 3u);
  for(const auto& m : matches)
    for(const auto& p : m.clippedPoly3D)
      EXPECT_TRUE(IsVertexOf(p, subject));
  for(const auto& m : inverted)
    for(const auto& p : m.clippedPoly3D)
      EXPECT_TRUE(IsVertexOf(p, clipped));
  EXPECT_GT(matcher.SnapRefused(), 0);
}

_FS_END_NAMESPACE
