#ifndef FSTOPOLOGYBUILDER_H
#define FSTOPOLOGYBUILDER_H

#include "FSClipping/FSFaceMatcher.h"
#include <FSMesh.h>
#include <FSMeshPolyFaceStorage.h>
#include <FSRegister.h>


_FS_BEGIN_NAMESPACE

struct FSTopologyData {
  std::vector<FSVec3> globalCoords;
  FSIntRegisterT cell2Node;
  FSMeshPolyFaceStorage polyFaces;
};

class FSTopologyBuilder
{

public:
  explicit FSTopologyBuilder(FS_floatT tol) : tol_(tol) {};

  FSTopologyData Build(const std::vector<FSFaceMatch>& matches);

private:
  FS_floatT tol_;
};


_FS_END_NAMESPACE

#endif