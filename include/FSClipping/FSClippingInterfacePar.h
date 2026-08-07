#ifndef FSCLIPPINGINTERFACEPAR_H
#define FSCLIPPINGINTERFACEPAR_H

#include <FSMesh.h>

#include "FSClipping/FSClippingEngine.h"


_FS_BEGIN_NAMESPACE

class FSClippingInterfacePar
{
public:
  FSClippingInterfacePar(FSClac& globalClac, FSClac& clac, FS_floatT tol, FS_intT boundaryMarkerMesh, FS_intT meshID)
    : globalClac_(globalClac), clac_(clac), tol_(tol), boundaryMarkerMesh_(boundaryMarkerMesh), meshID_(meshID){};

  FSMesh BuildSurfaceInterface(FSMesh& mesh, FS_intT meshID);

  FSMesh BuildVolumeInterface(FSMesh& mesh, FS_intT meshID);

private:
  FSMesh BuildInterface(FSMesh& mesh, FSClippingEngine::Mode mode, const FS_intT meshId);

  FSClac& globalClac_;
  FSClac& clac_;
  FS_floatT tol_;
  FS_intT meshID_;
  FS_intT boundaryMarkerMesh_;
};


_FS_END_NAMESPACE

#endif