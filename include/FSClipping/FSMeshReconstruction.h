#ifndef FSMESHRECONSTRUCTION_H
#define FSMESHRECONSTRUCTION_H

#include "FSClipping/FSTopologyData.h"
#include <FSMesh.h>
#include <FSUnstructMeshData.h>
#include <map>
#include <set>


_FS_BEGIN_NAMESPACE

// Shape of the clipped mesh: key = cell pool to create, value = its attribute names
using FSOutputSchema = std::map<FSMeshEnums::CellType, std::set<FSString>>;

class FSMeshReconstruction
{
public:
  explicit FSMeshReconstruction(FSClac& clac, FS_floatT tol = 1e-12) : clac_(clac), tol_(tol){};

  FSMesh Build(const FSUnstructMeshData& meshDataOriginal, FSTopologyData& topologyData);

  // Local shape only: for a single proc, or when the caller built every pool itself
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

  // Pools and attribute names seen on this proc, from the cells kept and their parents
  static FSOutputSchema LocalSchema(const FSUnstructMeshData& meshData, const FSTopologyData& topo);

  // Union of the schema over the procs: every proc then runs the same collective pool and attribute init
  void AgreeOnSchema(FSOutputSchema& schema);

  void CopyCellAttributes(const FSUnstructMeshData& meshData, const FSTopologyData& topo,
                          const FSOutputSchema& schema, FSUnstructMeshData& meshDataNew);

  FSClac& clac_;
  FS_floatT tol_;
};


_FS_END_NAMESPACE

#endif