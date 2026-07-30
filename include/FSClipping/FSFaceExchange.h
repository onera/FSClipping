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

// GatheredFaces aggregates all faces from one mesh group.
// cellToGlobalProc maps each volume cell to the global proc that owns it,
// allowing ScatterSend to route matches back without storing per-face proc IDs.
// localToGlobal maps sub-comm proc IDs to global proc IDs (populated even for
// procs that sent 0 faces so that ScatterSend can send empty match lists).
//
// The cell key is the (local proc ID, cell ID) pair, NOT the bare cell ID: cell
// IDs are local to a proc, so the several procs of one mesh may each use the same
// ID for a different cell. Keying on the ID alone let the last proc gathered
// overwrite the others, and every match for that ID was then scattered to a
// single wrong proc — leaving the other procs' cells unclipped and merging two
// distinct cells into one over-sized Poly3D downstream.
struct GatheredFaces {
  ReceivedFaces faces;
  // (owner's local proc ID, volume cell ID) → global proc ID
  std::map<std::pair<FS_intT, FS_intT>, FS_intT> cellToGlobalProc;
  std::map<FS_intT, FS_intT> localToGlobal; // local proc ID → global proc ID
};

// Clipper receives from both mesh groups in one call.
struct AllGatheredFaces {
  GatheredFaces meshA; // faces from procs that sent meshID == 1
  GatheredFaces meshB; // faces from procs that sent meshID == 2
};

// Pack face geometry + FSDM connectivity into clac's send buffer and transmit
// to destProc. destProc must call Receive(clac, thisProc).
void Send(FSClac& clac, FS_intT destProc, const std::vector<FSClippingFace>& faces);

void GatherSend(FSClac& clac, FS_intT clipperProc, FS_intT meshID, const std::vector<FSClippingFace>& faces);

// Receive and unpack faces and their FSDM connectivity from sourceProc.
ReceivedFaces Receive(FSClac& clac, FS_intT sourceProc);

AllGatheredFaces GatherReceiveAll(FSClac& clac);

} // namespace FSFaceExchange

namespace FSMatchExchange {

// Pack face geometry + FSDM connectivity into clac's send buffer and transmit
// to destProc. destProc must call Receive(clac, thisProc).
void Send(FSClac& clac, FS_intT destProc, const std::vector<FSFaceMatch>& matches);

// Receive and unpack faces and their FSDM connectivity from sourceProc.
std::vector<FSFaceMatch> Receive(FSClac& clac, FS_intT sourceProc);

static FSFaceMatch InvertSingle(const FSFaceMatch& m);

void ScatterSend(FSClac& clac, const std::vector<FSFaceMatch>& matches, const FSFaceExchange::GatheredFaces& gatheredA);

std::vector<FSFaceMatch> ScatterReceive(FSClac& clac, FS_intT clipperProc);

} // namespace FSMatchExchange

_FS_END_NAMESPACE

#endif // FSFACEEXCHANGE_H
