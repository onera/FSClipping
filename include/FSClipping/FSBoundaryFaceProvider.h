#ifndef FSBOUNDARYFACEPROVIDER_H
#define FSBOUNDARYFACEPROVIDER_H

#include <FSMesh.h>

#include "FSClipping/FSClippingFace.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

_FS_BEGIN_NAMESPACE

struct BoundaryExtraction {
  std::vector<FSClippingFace> faces;
  std::unordered_map<FS_intT, std::vector<FS_intT> > cell2BoundaryFaces; // cell ID -> local faces on the boundary
  std::unordered_set<FS_intT> surfaceCells;                               // local surface cells on the boundary (IDs unique across types)
};

class FSBoundaryFaceProvider
{
public:
  static BoundaryExtraction Extract(FSMesh& mesh, FSMeshFaceExtractor& ex, FS_intT boundaryMarker);
};


_FS_END_NAMESPACE

#endif