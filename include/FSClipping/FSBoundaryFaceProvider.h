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
};


_FS_END_NAMESPACE

#endif