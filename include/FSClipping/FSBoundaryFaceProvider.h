#ifndef FSBOUNDARYFACEPROVIDER_H
#define FSBOUNDARYFACEPROVIDER_H

#include <FSMesh.h>

#include "FSClipping/FSClippingFace.h"

_FS_BEGIN_NAMESPACE

class FSBoundaryFaceProvider
{
public:
  static std::vector<FSClippingFace> Extract(FSMesh& mesh,
                                             FSMeshFaceExtractor& ex,
                                             FS_intT boundaryMarker);

  static std::unordered_map<FS_intT, std::set<FS_intT> > ExtractOwnerCells(const std::vector<FSClippingFace>& faces);
};


_FS_END_NAMESPACE

#endif