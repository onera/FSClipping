#include "gtest/gtest.h"

#include "FSClipping/FSClippingUtil.h"

_FS_BEGIN_NAMESPACE

TEST(FSClippingTestGeo, AreFacesDifferent) {

  FSFloatArrayT coordinates1(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT coordArray1[FSMESH_NNODES_QUAD4][FS_3D] = {
      {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};

  FSFloatArrayT coordinates2(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT coordArray2[FSMESH_NNODES_QUAD4][FS_3D] = {
      {0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}};

  FSFloatArrayT coordinates3(FSMESH_NNODES_QUAD4, FS_3D);
  FS_floatT coordArray3[FSMESH_NNODES_QUAD4][FS_3D] = {
      {0, 0.5, 0}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}};

  for (FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i) {
    for (FS_intT d = 0; d < FS_3D; ++d) {
      coordinates1(i, d) = coordArray1[i][d];
      coordinates2(i, d) = coordArray2[i][d];
      coordinates3(i, d) = coordArray3[i][d];
    }
  }

  FSClippingFace f1(coordinates1);
  FSClippingFace f2(coordinates2);
  FSClippingFace f3(coordinates3);

  EXPECT_TRUE(FSClippingUtil::AreFacesIdentical(f1, f2, 1e-12));
  EXPECT_FALSE(FSClippingUtil::AreFacesIdentical(f2, f3, 1e-12));
}

// ---- 1. Identical squares (currently commented out)
TEST(FSClippingTestGeo, ClipConvexPolygon_IdenticalSquares) {
  {
    std::vector<FSVec2> square = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    std::vector<FSVec2> out;
    EXPECT_TRUE(FSClippingUtil::ClipConvexPolygon(square, square, out));
    EXPECT_EQ(4, out.size());
  }
}

// ---- 2. Full inclusion
TEST(FSClippingTestGeo, ClipConvexPolygon_FullInclusion) {
  std::vector<FSVec2> big = {{-1, -1}, {2, -1}, {2, 2}, {-1, 2}};
  std::vector<FSVec2> small = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  std::vector<FSVec2> out;
  EXPECT_TRUE(FSClippingUtil::ClipConvexPolygon(small, big, out));
  EXPECT_TRUE(out.size() == 4);
}

// ---- 3. Partial overlap (classic test)
TEST(FSClippingTestGeo, ClipConvexPolygon_PartialOverlap) {
  std::vector<FSVec2> P = {{0, 0}, {2, 0}, {2, 1}, {0, 1}};
  std::vector<FSVec2> Q = {{1, -1}, {3, -1}, {3, 2}, {1, 2}};
  std::vector<FSVec2> out;
  std::vector<FSVec2> result = {{1, 0}, {2, 0}, {2, 1}, {1, 1}};

  bool ok = FSClippingUtil::ClipConvexPolygon(P, Q, out);
  ASSERT_TRUE(ok);
  // Expected shape : rectangle [1,0]x[2,1] → 4 pts
  EXPECT_TRUE(out.size() == 4);
  for (std::size_t i = 0; i < result.size(); ++i)
    ASSERT_TRUE((out[i] - result[i]).Norm() < 1e-12);
}

// ---- 4. No overlap
TEST(FSClippingTestGeo, ClipConvexPolygon_NoOverlap) {
  std::vector<FSVec2> A = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  std::vector<FSVec2> B = {{2, 0}, {3, 0}, {3, 1}, {2, 1}};
  std::vector<FSVec2> out;

  bool ok = FSClippingUtil::ClipConvexPolygon(A, B, out);
  ASSERT_FALSE(ok);
  ASSERT_TRUE(out.empty());
}

// ---- 5. Triangle clipped against square
TEST(FSClippingTestGeo, ClipConvexPolygon_TriangleVsSquare) {
  std::vector<FSVec2> tri = {{0, 0}, {2, 0}, {1, 2}};
  std::vector<FSVec2> sq = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  std::vector<FSVec2> out;
  std::vector<FSVec2> result = {{0, 0}, {1, 0}, {1, 1}, {0.5, 1}};

  bool ok = FSClippingUtil::ClipConvexPolygon(tri, sq, out);
  ASSERT_TRUE(ok);
  ASSERT_TRUE(out.size() >= 3);

  for (std::size_t i = 0; i < result.size(); ++i)
    ASSERT_TRUE((out[i] - result[i]).Norm() < 1e-12);
}

_FS_END_NAMESPACE
