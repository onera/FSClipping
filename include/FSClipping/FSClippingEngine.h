#ifndef FSCLIPPINGENGINE_H
#define FSCLIPPINGENGINE_H

#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSTopologyAssembler.h"
#include <FSMesh.h>

_FS_BEGIN_NAMESPACE

// Shared clipping pipeline used by every public entry point.
class FSClippingEngine
{
public:
  enum class Mode { Surface, Volume };

  explicit FSClippingEngine(FS_floatT tol) : tol_(tol) {}

  FSTopologyData BuildTopology(FSMesh& mesh, const BoundaryExtraction& be, const std::vector<FSFaceMatch>& matches,
                               Mode mode);

  FSMesh Reconstruct(FSClac& clac, FSMesh& meshOriginal, FSTopologyData& topo, bool copyAttributes = false);

  static bool RunMatcherProc(FSClac& globalClac, FSClac& localClac, FS_floatT tol);

private:
  FS_floatT tol_;
};

_FS_END_NAMESPACE

#endif
