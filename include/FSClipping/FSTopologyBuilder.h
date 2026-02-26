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
  FSIntRegisterT cell2Node;
  FSMeshPolyFaceStorage polyFaces;
};

class FSTopologyBuilder
{

public:
  explicit FSTopologyBuilder(FS_floatT tol) : tol_(tol),
                                              cell2NodeBuilder_(std::make_unique<FSCell2NodeBuilder>(FSCell2NodeBuilder(tol))) {};

  FSTopologyData BuildSurfaceTopo(const std::vector<FSFaceMatch>& matches);
  FSTopologyData BuildVolumeTopo(const FSIntArrayT& cell2Node,
                                 const std::set<FS_intT>& bdryCells,
                                 const FSCellPool& cellPool,
                                 const FSFloatArrayT& oldCoords);
  //  FSTopologyData Add

private:
  FS_floatT tol_;

  std::unique_ptr<FSCell2NodeBuilder> cell2NodeBuilder_;         // On souhaiterai se passer de la variable cell2NodeBuilder_ :
                                                                 //  0- ajouter tout les precedents noeuds du matches, donc faire un constructeur avec TopologyData
                                                                 //  1- Dans AddVolumeCellNodes -> garder seulement le celldData en memoire +
  std::unique_ptr<FSPolyFaceBuilder> polyFaceBuilder_ = nullptr; // On souhaiterai aussi se passer de polyFaceBuilder
};


_FS_END_NAMESPACE

#endif