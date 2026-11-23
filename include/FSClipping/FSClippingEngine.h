#ifndef FSCLIPPINGENGINE_H
#define FSCLIPPINGENGINE_H

#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSNodeSnap.h"
#include "FSClipping/FSTopologyData.h"
#include <FSMesh.h>

_FS_BEGIN_NAMESPACE

// Shared clipping pipeline used by every public entry point.
class FSClippingEngine
{
public:
  explicit FSClippingEngine(FS_floatT tol) : tol_(tol) {}

  // `moves` relocates old nodes, owned and ghost, before any vertex is resolved
  FSTopologyData BuildTopology(FSMesh& mesh, const BoundaryExtraction& be, const std::vector<FSFaceMatch>& matches,
                               const std::vector<FSNodeMove>& moves = {});

  FSMesh Reconstruct(FSClac& clac, FSMesh& meshOriginal, FSTopologyData& topo, bool copyAttributes = false);

  // snapTol > 0: mesh 2 interface nodes within snapTol of a mesh 1 node are moved onto it
  static bool RunMatcherProc(FSClac& globalClac, FSClac& localClac, FS_floatT tol, FS_floatT snapTol = 0.0);

private:
  FS_floatT tol_;
};

_FS_END_NAMESPACE

#endif
