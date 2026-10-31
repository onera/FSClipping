#ifndef FSTOPOLOGYDATA_H
#define FSTOPOLOGYDATA_H

#include <FSMesh.h>
#include <FSMeshPolyFaceStorage.h>
#include <FSRegister.h>
#include <FSVec3.h>
#include <vector>

_FS_BEGIN_NAMESPACE

// Contract between the topology reconstruction and FSMeshReconstruction
struct FSTopologyData {
  std::vector<FSVec3> globalCoords;
  FSIntRegisterT cell2NodePoly3D;
  FSIntRegisterT cell2NodePoly2D;
  FSCellType2IntArrayT cell2NodeInner;
  FSMeshPolyFaceStorage polyFaces;

  FSCellType2IntArrayT cellParent;
  FSCellType2IntArrayT cellParentType;

  // Global numbering assigned by the clipper proc
  FSIntArrayT nodeGlobalNumbers;   // parallel to globalCoords

  // CADGroupID of each Poly2D
  FSIntArrayT poly2DMarkers;
};

_FS_END_NAMESPACE

#endif
