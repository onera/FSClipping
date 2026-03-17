#ifndef FSTOPOLOGYBUILDER_H
#define FSTOPOLOGYBUILDER_H

#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSPolyFaceBuilder.h"
#include <FSMesh.h>
#include <FSMeshPolyFaceStorage.h>
#include <FSRegister.h>
#include <memory>


_FS_BEGIN_NAMESPACE

struct FSTopologyData {
  std::vector<FSVec3> globalCoords;
  FSIntRegisterT cell2NodePoly3D;
  FSIntRegisterT cell2NodePoly2D;
  FSCellType2IntArrayT cell2NodeInner;
  FSMeshPolyFaceStorage polyFaces;
};

class FSTopologyAssembler
{

public:
  explicit FSTopologyAssembler(FS_floatT tol) : tol_(tol),
                                                cell2NodeBuilder_(FSCell2NodeBuilder(tol)) {};

  FSTopologyData BuildSurfaceTopo(const std::vector<FSFaceMatch>& matches);

  void BuildVolumeTopo(const FSIntArrayT& cell2Node,
                       const std::set<FS_intT>& bdryCells,
                       const FSCellPool& cellPool,
                       const FSFloatArrayT& oldCoords,
                       FSTopologyData& topology);

  FSIntArrayT UpdateOldCell2Node(const FSIntArrayT& oldCell2Node,
                                 const FSFloatArrayT& oldCoords);

private:
  void CheckSurfaceWasBuilt() const;

  FS_floatT tol_;
  FSCell2NodeBuilder cell2NodeBuilder_;                          // On souhaiterai se passer de la variable cell2NodeBuilder_ :
                                                                 //  0- ajouter tout les precedents noeuds du matches, donc faire un constructeur avec TopologyData
                                                                 //  1- Dans AddVolumeCellNodes -> garder seulement le celldData en memoire +
  std::unique_ptr<FSPolyFaceBuilder> polyFaceBuilder_ = nullptr; // On souhaiterai aussi se passer de polyFaceBuilder et fonctionner seulement avec un FSTopologyData(Surface)
  bool surfaceBuilt_ = false;
};


_FS_END_NAMESPACE

#endif