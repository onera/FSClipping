#ifndef FSTOPOLOGYASSEMBLER_H
#define FSTOPOLOGYASSEMBLER_H

#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSPolyFaceBuilder.h"
#include <FSMesh.h>
#include <FSMeshPolyFaceStorage.h>
#include <FSRegister.h>
#include <memory>


_FS_BEGIN_NAMESPACE

struct CellParentInfo {
  FSMeshEnums::CellType parentType;
  FS_intT parentId;
};

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
  FSIntArrayT poly2DGlobalNumbers; // parallel to the Poly2D cells (match order)

  // CADGroupID of each Poly2D
  FSIntArrayT poly2DMarkers;

  void Send(FSClac& clac, FS_intT destProc);
  void Received(FSClac& clac, FS_intT sourceProc);
};



class FSTopologyAssembler
{

public:
  explicit FSTopologyAssembler(FS_floatT tol) : tol_(tol), cell2NodeBuilder_(FSCell2NodeBuilder(tol)){};

  FSTopologyData BuildSurfaceTopo(const std::vector<FSFaceMatch>& matches,
                                  const std::unordered_set<GeomFaceKey, GeomFaceKeyHash>& faceKeys);

  void BuildVolumeTopo(const FSIntArrayT& cell2Node, const std::set<FS_intT>& bdry3DCells, const FSCellPool& cellPool,
                       const FSFloatArrayT& oldCoords, FSTopologyData& topology);

  void AppendUnclippedSurfaces(FSMesh& mesh, const std::unordered_map<FS_intT, std::set<FS_intT> >& bdry2DCells,
                               const FSFloatArrayT& oldCoords, FSTopologyData& topo);

  FSIntArrayT UpdateOldCell2Node(const FSIntArrayT& oldCell2Node, const FSCellPool& cellPool,
                                 const FSFloatArrayT& oldCoords, const std::set<FS_intT>& bdry2DCells,
                                 FSIntArrayT& cellParent);

private:
  void CheckSurfaceWasBuilt() const;

  FS_floatT tol_;
  FSCell2NodeBuilder cell2NodeBuilder_;
  std::unique_ptr<FSPolyFaceBuilder> polyFaceBuilder_ = nullptr;
  bool surfaceBuilt_ = false;
};


_FS_END_NAMESPACE

#endif