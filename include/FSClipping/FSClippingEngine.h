#ifndef FSCLIPPINGENGINE_H
#define FSCLIPPINGENGINE_H

#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSTopologyAssembler.h"
#include <FSMesh.h>

_FS_BEGIN_NAMESPACE

//! Shared clipping pipeline used by every public entry point
//! (FSClippingInterface, FSClippingInterfacePar, FSClippedMesh).
//!
//! The engine owns the steps that are independent of where the meshes and
//! the matches come from: topology assembly and mesh reconstruction.
//! Boundary extraction and match acquisition (local matcher vs MPI
//! exchange) stay in the calling interfaces.
class FSClippingEngine
{
public:
  enum class Mode { Surface, Volume };

  explicit FSClippingEngine(FS_floatT tol) : tol_(tol) {}

  //! Build the clipped topology of mesh from its boundary extraction and
  //! the pre-computed face matches.
  //! Surface mode: clipped boundary faces only.
  //! Volume mode: clipped boundary faces + volume cells + unclipped surfaces.
  FSTopologyData BuildTopology(FSMesh& mesh, const BoundaryExtraction& be, const std::vector<FSFaceMatch>& matches,
                               Mode mode);

  //! Reconstruct the final mesh from the topology data.
  //! copyAttributes also transfers the original mesh attributes (used by
  //! the DataManagerOp entry point).
  FSMesh Reconstruct(FSClac& clac, FSMesh& meshOriginal, FSTopologyData& topo, bool copyAttributes = false);

  //! Role of the matcher proc (proc 2) in the 3-proc parallel layout:
  //! receive the boundary faces of both meshes, compute the matches and
  //! send them back to each owner proc.
  static void RunMatcherProc(FSClac& globalClac, FSClac& localClac, FS_floatT tol);

private:
  FS_floatT tol_;
};

_FS_END_NAMESPACE

#endif
