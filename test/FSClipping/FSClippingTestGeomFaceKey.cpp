#include "gtest/gtest.h"

#include "FSClipping/FSBoundaryFaceProvider.h"

#include <unordered_set>

_FS_BEGIN_NAMESPACE

TEST(FSClippingTestGeomFaceKey, QuantizedPointEquality)
{
  QuantizedPoint p1{1, 2, 3};
  QuantizedPoint p2{1, 2, 3};
  QuantizedPoint p3{1, 2, 4};

  EXPECT_TRUE(p1 == p2);
  EXPECT_FALSE(p1 == p3);
}

TEST(FSClippingTestGeomFaceKey, QuantizedPointLessThan)
{
  // Test lexicographic ordering used for canonicalization
  QuantizedPoint p1{1, 2, 3};
  QuantizedPoint p2{1, 2, 4};
  QuantizedPoint p3{1, 3, 0};
  QuantizedPoint p4{2, 0, 0};

  EXPECT_TRUE(p1 < p2); // same x, same y, z1 < z2
  EXPECT_TRUE(p2 < p3); // same x, y1 < y2
  EXPECT_TRUE(p3 < p4); // x1 < x2
}

TEST(FSClippingTestGeomFaceKey, CanonicalizeRotationInvariance)
{
  // Triangle: test rotation invariance
  FSFloatArrayT coords1(3, FS_3D);
  coords1(0, 0) = 0.0;
  coords1(0, 1) = 0.0;
  coords1(0, 2) = 0.0;
  coords1(1, 0) = 1.0;
  coords1(1, 1) = 0.0;
  coords1(1, 2) = 0.0;
  coords1(2, 0) = 0.0;
  coords1(2, 1) = 1.0;
  coords1(2, 2) = 0.0;

  // Rotated: [1,0,0] -> [0,1,0] -> [0,0,0]
  FSFloatArrayT coords2(3, FS_3D);
  coords2(0, 0) = 1.0;
  coords2(0, 1) = 0.0;
  coords2(0, 2) = 0.0;
  coords2(1, 0) = 0.0;
  coords2(1, 1) = 1.0;
  coords2(1, 2) = 0.0;
  coords2(2, 0) = 0.0;
  coords2(2, 1) = 0.0;
  coords2(2, 2) = 0.0;

  GeomFaceKey key1(coords1, 1.0);
  GeomFaceKey key2(coords2, 1.0);

  EXPECT_TRUE(key1 == key2);
}

TEST(FSClippingTestGeomFaceKey, CanonicalizeReversalInvariance)
{
  // Test that reversed order (opposite winding) gives the same canonical key
  FSFloatArrayT coords1(3, FS_3D);
  coords1(0, 0) = 0.0;
  coords1(0, 1) = 0.0;
  coords1(0, 2) = 0.0;
  coords1(1, 0) = 1.0;
  coords1(1, 1) = 0.0;
  coords1(1, 2) = 0.0;
  coords1(2, 0) = 0.0;
  coords1(2, 1) = 1.0;
  coords1(2, 2) = 0.0;

  // Reversed: [0,1,0] -> [1,0,0] -> [0,0,0]
  FSFloatArrayT coords2(3, FS_3D);
  coords2(0, 0) = 0.0;
  coords2(0, 1) = 1.0;
  coords2(0, 2) = 0.0;
  coords2(1, 0) = 1.0;
  coords2(1, 1) = 0.0;
  coords2(1, 2) = 0.0;
  coords2(2, 0) = 0.0;
  coords2(2, 1) = 0.0;
  coords2(2, 2) = 0.0;

  GeomFaceKey key1(coords1, 1.0);
  GeomFaceKey key2(coords2, 1.0);

  EXPECT_TRUE(key1 == key2);
}

TEST(FSClippingTestGeomFaceKey, CanonicalizeQuadWithReversal)
{
  // Square: test with 4 vertices and reversal
  FSFloatArrayT coords1(4, FS_3D);
  coords1(0, 0) = 0.0;
  coords1(0, 1) = 0.0;
  coords1(0, 2) = 0.0;
  coords1(1, 0) = 1.0;
  coords1(1, 1) = 0.0;
  coords1(1, 2) = 0.0;
  coords1(2, 0) = 1.0;
  coords1(2, 1) = 1.0;
  coords1(2, 2) = 0.0;
  coords1(3, 0) = 0.0;
  coords1(3, 1) = 1.0;
  coords1(3, 2) = 0.0;

  // Reversed order (counter-clockwise instead of clockwise)
  FSFloatArrayT coords2(4, FS_3D);
  coords2(0, 0) = 0.0;
  coords2(0, 1) = 1.0;
  coords2(0, 2) = 0.0;
  coords2(1, 0) = 1.0;
  coords2(1, 1) = 1.0;
  coords2(1, 2) = 0.0;
  coords2(2, 0) = 1.0;
  coords2(2, 1) = 0.0;
  coords2(2, 2) = 0.0;
  coords2(3, 0) = 0.0;
  coords2(3, 1) = 0.0;
  coords2(3, 2) = 0.0;

  GeomFaceKey key1(coords1, 1.0);
  GeomFaceKey key2(coords2, 1.0);

  EXPECT_TRUE(key1 == key2);
}

TEST(FSClippingTestGeomFaceKey, HashSameKeySameHash)
{
  // Equal keys must have equal hashes (required for unordered_set correctness)
  FSFloatArrayT coords1(3, FS_3D);
  coords1(0, 0) = 0.0;
  coords1(0, 1) = 0.0;
  coords1(0, 2) = 0.0;
  coords1(1, 0) = 1.0;
  coords1(1, 1) = 0.0;
  coords1(1, 2) = 0.0;
  coords1(2, 0) = 0.0;
  coords1(2, 1) = 1.0;
  coords1(2, 2) = 0.0;

  FSFloatArrayT coords2(3, FS_3D);
  coords2(0, 0) = 1.0;
  coords2(0, 1) = 0.0;
  coords2(0, 2) = 0.0;
  coords2(1, 0) = 0.0;
  coords2(1, 1) = 1.0;
  coords2(1, 2) = 0.0;
  coords2(2, 0) = 0.0;
  coords2(2, 1) = 0.0;
  coords2(2, 2) = 0.0;

  GeomFaceKey key1(coords1, 1.0);
  GeomFaceKey key2(coords2, 1.0);

  GeomFaceKeyHash hasher;

  EXPECT_TRUE(key1 == key2);
  EXPECT_EQ(hasher(key1), hasher(key2));
}

TEST(FSClippingTestGeomFaceKey, UnorderedSetRotationInvariance)
{
  // Verify that rotated faces are correctly identified as duplicates
  std::unordered_set<GeomFaceKey, GeomFaceKeyHash> faceKeys;

  FSFloatArrayT coords1(3, FS_3D);
  coords1(0, 0) = 0.0;
  coords1(0, 1) = 0.0;
  coords1(0, 2) = 0.0;
  coords1(1, 0) = 1.0;
  coords1(1, 1) = 0.0;
  coords1(1, 2) = 0.0;
  coords1(2, 0) = 0.0;
  coords1(2, 1) = 1.0;
  coords1(2, 2) = 0.0;

  FSFloatArrayT coords2(3, FS_3D);
  coords2(0, 0) = 1.0;
  coords2(0, 1) = 0.0;
  coords2(0, 2) = 0.0;
  coords2(1, 0) = 0.0;
  coords2(1, 1) = 1.0;
  coords2(1, 2) = 0.0;
  coords2(2, 0) = 0.0;
  coords2(2, 1) = 0.0;
  coords2(2, 2) = 0.0;

  faceKeys.insert(GeomFaceKey(coords1, 1.0));

  // key2 should be found in the set since it's equal to key1 (rotationally equivalent)
  EXPECT_TRUE(faceKeys.find(GeomFaceKey(coords2, 1.0)) != faceKeys.end());
  EXPECT_EQ(faceKeys.size(), 1); // Only one unique face
}

_FS_END_NAMESPACE