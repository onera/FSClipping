#ifndef FSFACEEXCHANGE_H
#define FSFACEEXCHANGE_H

#include "FSClipping/FSClippingFace.h"
#include "FSClipping/FSFaceMatcher.h"
#include <FSClac.h>
#include <FSHashableFace.h>
#include <vector>

_FS_BEGIN_NAMESPACE

// Exchange of FSClippingFace geometry and FSDM connectivity via FSClac
namespace FSFaceExchange {

struct ReceivedFaces {
  std::vector<FSClippingFace> faces;
  std::vector<FSFaceConnectivity> connectivity;
};

struct GatheredFaces {
  ReceivedFaces faces;
  std::map<std::pair<FS_intT, FS_intT>, FS_intT> cellToGlobalProc; // <local proc ID, volume cell ID> → global proc ID
  std::map<FS_intT, FS_intT> localToGlobal;                        // local proc ID → global proc ID
};

struct AllGatheredFaces {
  GatheredFaces meshA;
  GatheredFaces meshB;
};

// Send all the faces (faceIndex, coord[3D], connectivity) of each proc to the only clipperProc (=0)
void GatherSend(FSClac& clac, FS_intT clipperProc, FS_intT meshID, const std::vector<FSClippingFace>& faces);

// Receive all the faces of each proc (GatherSend) and stored them in the AllGatheredFaces struct
AllGatheredFaces GatherReceiveAll(FSClac& clac);

} // namespace FSFaceExchange

namespace FSMatchExchange {

void Send(FSClac& clac, FS_intT destProc, const std::vector<FSFaceMatch>& matches);

void ScatterSend(FSClac& clac, const std::vector<FSFaceMatch>& matches, const FSFaceExchange::GatheredFaces& gatheredA);

std::vector<FSFaceMatch> ScatterReceive(FSClac& clac, FS_intT clipperProc);

} // namespace FSMatchExchange

_FS_END_NAMESPACE

#endif // FSFACEEXCHANGE_H
