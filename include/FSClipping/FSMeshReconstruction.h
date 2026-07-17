#ifndef FSMESHRECONSTRUCTION_H
#define FSMESHRECONSTRUCTION_H

#include "FSClipping/FSTopologyAssembler.h"
#include <FSMesh.h>


_FS_BEGIN_NAMESPACE

class FSMeshReconstruction
{
public:
  explicit FSMeshReconstruction(FSClac& clac) : clac_(clac){};

  FSMesh Build(const FSUnstructMeshData& meshDataOriginal, FSTopologyData& topologyData);

  void CopyCellAttributes(const FSUnstructMeshData& meshData, const FSTopologyData& topo,
                          FSUnstructMeshData& meshDataNew);

  void CopyAttributes(const FSMesh& meshOriginal, FSMesh& meshClipped);

private:
  //! Level 2 parallel assembly of the distributed node numbering.
  /*!
    FSDM's global-numbering contract: each node is owned by exactly one proc
    (contiguous ranges per proc) and cell2node references those global numbers.
    This routine uses the clipper-assigned IDs (topo.nodeGlobalNumbers) as merge
    key to (1) decide ownership of nodes shared across procs, (2) compute the
    contiguous distributed numbering, (3) remap cell2NodePoly2D/Poly3D in place.
    \param[in,out] topo             topology whose connectivity gets remapped.
    \param[out] ownedLocalRows      local row indices (into globalCoords) of the owned nodes, in contiguous-ID order.
    \param[out] ownedGlobalNumbers  clipper IDs of the owned nodes (stable GlobalNumber attribute).
  */
  void RemapToDistributedNumbering(FSTopologyData& topo, std::vector<FS_intT>& ownedLocalRows,
                                   FSIntArrayT& ownedGlobalNumbers);

  FSClac& clac_;
};


_FS_END_NAMESPACE

#endif