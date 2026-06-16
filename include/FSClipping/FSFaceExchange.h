#ifndef FSFACEEXCHANGE_H
#define FSFACEEXCHANGE_H

#include "FSClipping/FSClippingFace.h"
#include "FSClipping/FSFaceMatcher.h"
#include <FSClac.h>
#include <FSHashableFace.h>
#include <vector>

_FS_BEGIN_NAMESPACE

// Point-to-point exchange of FSClippingFace geometry and FSDM connectivity
// via the FSClac buffer protocol.
//
// Send transmits both vertex geometry and the per-face FSFaceConnectivity
// (owner cell + neighbor face IDs and types). The received connectivity is
// used to populate FSFaceMatch topology fields after ComputeMatches so that
// BuildSurfaceTopo can be called on the clipper proc.
//
// Precondition for Send: faces must come from FSBoundaryFaceProvider::Extract
// (valid _faceFSDM). Geometry-only faces (null topology) are not supported.
namespace FSFaceExchange {

struct ReceivedFaces {
  std::vector<FSClippingFace> faces;
  std::vector<FSFaceConnectivity> connectivity; // parallel to faces
};

// Pack face geometry + FSDM connectivity into clac's send buffer and transmit
// to destProc. destProc must call Receive(clac, thisProc).
void Send(FSClac& clac, FS_intT destProc, const std::vector<FSClippingFace>& faces);

// Receive and unpack faces and their FSDM connectivity from sourceProc.
ReceivedFaces Receive(FSClac& clac, FS_intT sourceProc);

} // namespace FSFaceExchange

namespace FSMatchExchange {

// Pack face geometry + FSDM connectivity into clac's send buffer and transmit
// to destProc. destProc must call Receive(clac, thisProc).
void Send(FSClac& clac, FS_intT destProc, const std::vector<FSFaceMatch>& matches);

// Receive and unpack faces and their FSDM connectivity from sourceProc.
std::vector<FSFaceMatch> Receive(FSClac& clac, FS_intT sourceProc);

} // namespace FSMatchExchange

_FS_END_NAMESPACE

#endif // FSFACEEXCHANGE_H
