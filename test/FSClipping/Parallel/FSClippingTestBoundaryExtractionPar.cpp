// Partition-independence of the boundary extraction, run at exactly 4 MPI procs.
//
// The number of faces extracted on a boundary marker is a property of the mesh, so summing it over all procs must
// give the same total whatever the partitioning. RCB partitions volume and surface cells independently, so a
// hexahedron and the quad closing it can land on different procs. Such a pair is reported twice by the face extractor
// -- once on the volume's proc (surface neighbour remote) and once on the surface's proc (surface cell as owner) --
// and must end up in the extraction exactly once.
//
// cube_hexa_coarse_par.grid has 16 faces on marker 6, and 4 procs is the smallest count at which RCB actually splits
// such a pair on this mesh.

#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceSeparator.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"

_FS_BEGIN_NAMESPACE

namespace {

/// Number of faces on marker 6 of cube_hexa_coarse_par.grid, as extracted on a single proc.
constexpr FS_int32T sExpectedFaceCount = 16;
constexpr FS_intT sMarker = 6;

FSMesh LoadRepartitionedCube(FSClac& clac)
{
  FSMesh mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse_par.grid"));
  polyMeshRepartition(&clac, mesh.GetMeshData());
  mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  return mesh;
}

} // namespace

TEST(FSClippingTestBoundaryExtractionPar, FaceCountIsPartitionIndependent)
{
  FSClac clac;
  SKIP_UNLESS_EXACT_PROCS(clac, 4);

  FSMesh mesh = LoadRepartitionedCube(clac);

  FSMeshFaceExtractor ex;
  const BoundaryExtraction extraction = FSBoundaryFaceProvider::Extract(mesh, ex, sMarker, 1e-8);

  FS_int32T localFaces = static_cast<FS_int32T>(extraction.faces.size());
  FS_int32T totalFaces = 0;
  clac.Sum(&localFaces, &totalFaces, 1);

  // A volume/surface pair split across two procs must still be extracted, otherwise the total silently drops below
  // the sequential reference.
  EXPECT_EQ(totalFaces, sExpectedFaceCount);
}

TEST(FSClippingTestBoundaryExtractionPar, NoFaceIsExtractedTwice)
{
  FSClac clac;
  SKIP_UNLESS_EXACT_PROCS(clac, 4);

  FSMesh mesh = LoadRepartitionedCube(clac);

  FSMeshFaceExtractor ex;
  const BoundaryExtraction extraction = FSBoundaryFaceProvider::Extract(mesh, ex, sMarker, 1e-8);

  // Each proc must hold distinct faces ...
  EXPECT_EQ(extraction.faceKeys.size(), extraction.faces.size()) << "a face was extracted twice on this proc";

  // ... and a split pair, reported to both of its procs, must not be kept by both.
  FS_int32T localKeys = static_cast<FS_int32T>(extraction.faceKeys.size());
  FS_int32T totalKeys = 0;
  clac.Sum(&localKeys, &totalKeys, 1);
  EXPECT_EQ(totalKeys, sExpectedFaceCount);
}

_FS_END_NAMESPACE
