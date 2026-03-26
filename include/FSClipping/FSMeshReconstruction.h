#ifndef FSMESHRECONSTRUCTION_H
#define FSMESHRECONSTRUCTION_H

#include "FSClipping/FSTopologyAssembler.h"
#include <FSMesh.h>


_FS_BEGIN_NAMESPACE

class FSMeshReconstruction
{
public:
  explicit FSMeshReconstruction(FSClac& clac) : clac_(clac) {};

  FSMeshData Build(const FSUnstructMeshData& meshDataOriginal, FSTopologyData& topologyData);

  void CopyCellAttributes(const FSUnstructMeshData& meshData,
                          const FSTopologyData& topo,
                          FSUnstructMeshData& meshDataNew);

private:
  FSClac& clac_;
};


_FS_END_NAMESPACE

#endif