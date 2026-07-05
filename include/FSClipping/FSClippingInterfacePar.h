#ifndef FSCLIPPINGINTERFACEPAR_H
#define FSCLIPPINGINTERFACEPAR_H

#include <FSMesh.h>

#include "FSClipping/FSClippingEngine.h"


_FS_BEGIN_NAMESPACE

class FSClippingInterfacePar
{
public:
  FSClippingInterfacePar(FSClac& globalClac, FSClac& clac, FS_floatT tol, FS_intT boundaryMarkerMesh)
    : globalClac_(globalClac), clac_(clac), tol_(tol), boundaryMarkerMesh_(boundaryMarkerMesh) {};

  FSMesh BuildSurfaceInterface(FSMesh& mesh);

  FSMesh BuildVolumeInterface(FSMesh& mesh);

private:
  FSMesh BuildInterface(FSMesh& mesh, FSClippingEngine::Mode mode, bool matchRemoteFaces);

  FSClac& globalClac_;
  FSClac& clac_;
  FS_floatT tol_;
  FS_intT boundaryMarkerMesh_;
};


_FS_END_NAMESPACE

#endif