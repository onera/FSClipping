#ifndef FSMESHRECONSTRUCTION_H
#define FSMESHRECONSTRUCTION_H

#include "FSClipping/FSTopologyAssembler.h"
#include <FSMesh.h>


_FS_BEGIN_NAMESPACE

class FSMeshReconstruction
{
public:
  explicit FSMeshReconstruction(FSClac& clac, FS_floatT tol = 1e-12) : clac_(clac), tol_(tol){};

  FSMesh Build(const FSUnstructMeshData& meshDataOriginal, FSTopologyData& topologyData);

  void CopyCellAttributes(const FSUnstructMeshData& meshData, const FSTopologyData& topo,
                          FSUnstructMeshData& meshDataNew);

  void CopyAttributes(const FSMesh& meshOriginal, FSMesh& meshClipped);

private:
  /*
    Parallel assembly of the distributed node numbering :
      1 - Decide ownership of nodes shared across procs
      2 - Compute the contiguous distributed numbering
      3 - Remap cell2NodePoly2D/Poly3D/cell2NodeInner in place.
    We do the remapping because each node is owned by exactly one proc (contiguous ranges per proc) and cell2node
    references those global numbers.
  */
  void RemapToDistributedNumbering(FSTopologyData& topo, std::vector<FS_intT>& ownedLocalRows,
                                   FSIntArrayT& ownedGlobalNumbers);

  FSClac& clac_;
  FS_floatT tol_;
};


_FS_END_NAMESPACE

#endif