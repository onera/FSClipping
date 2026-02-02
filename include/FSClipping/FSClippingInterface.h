#ifndef FSCLIPPINGINTERFACE_H
#define FSCLIPPINGINTERFACE_H

#include <FSMesh.h>

#include "FSClipping/FSClippingFace.h"

_FS_BEGIN_NAMESPACE

class FSClippingInterface
{
public:
  FSClippingInterface(FSClac& clac1,
                      FSClac& clac2,
                      FS_floatT tol,
                      FS_intT markerMesh1,
                      FS_intT markerMesh2) : clac1_(clac1), clac2_(clac2), tol_(tol), boundaryMarkerMesh1_(markerMesh1), boundaryMarkerMesh2_(markerMesh2) {};

  FSMeshData BuildInterface(FSMesh& mesh1, FSMesh& mesh2, FSClac& clacInterface);

private:
  FS_floatT tol_;
  FSClac& clac1_;
  FSClac& clac2_;
  FS_intT boundaryMarkerMesh1_;
  FS_intT boundaryMarkerMesh2_;
};


_FS_END_NAMESPACE

#endif