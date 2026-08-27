#include "FSClipping/FSFaceExchange.h"
#include <unordered_map>

_FS_BEGIN_NAMESPACE

namespace FSFaceExchange {

void GatherSend(FSClac& globalClac, FS_intT clipperProc, FS_intT meshID, const std::vector<FSClippingFace>& faces)
{
  FS_intT worldProcID = globalClac.GetWorldProcID();
  FS_intT localProcID = globalClac.GetProcID();

  // Header: worldProcID + localProcID + meshID + n_faces
  FSClac::sizeT bufSize = globalClac.GetBufSizeInt32(4);
  for(const auto& f : faces) {
    bufSize += globalClac.GetBufSizeInt32(2);
    bufSize += globalClac.GetBufSizeFloat64(static_cast<FSClac::intT>(f.vertices().size()) * 3);
    bufSize += f.topo()._faceFSDM->GetBufSize(globalClac);
  }

  globalClac.SelectSendBuffer(clipperProc);
  globalClac.InitSendBuffer(bufSize);
  globalClac.Pack(&worldProcID, 1);
  globalClac.Pack(&localProcID, 1);
  globalClac.Pack(&meshID, 1);
  FS_intT n = static_cast<FS_intT>(faces.size());
  globalClac.Pack(&n, 1);

  for(const auto& f : faces) {
    FS_intT nv = static_cast<FS_intT>(f.vertices().size());
    globalClac.Pack(&nv, 1);
    FS_intT faceIndex = f.faceIndex();
    globalClac.Pack(&faceIndex, 1);
    for(const auto& v : f.vertices()) {
      FS_float64T coords[3] = {v[0], v[1], v[2]};
      globalClac.Pack(coords, 3);
    }
    FSCellFace owner = f.topo()._faceFSDM->mOwner;
    FSCellFace neighbor = f.topo()._faceFSDM->mNeighbor;
    owner.Pack(globalClac);
    neighbor.Pack(globalClac);
  }

  globalClac.SendBuffer(clipperProc);
}

AllGatheredFaces GatherReceiveAll(FSClac& globalClac)
{
  // 1 - Receive raw data from every non-clipper proc.
  struct RawFace {
    FS_intT faceIndex;
    FSFloatArrayT coords;
    FSFaceConnectivity conn;
    FS_intT originGlobalProc;
    FS_intT meshID; // 1 = mesh A, 2 = mesh B
  };

  struct RawMessage {
    FS_intT worldProcID, localProcID, meshID;
    std::vector<RawFace> faces;
  };

  const FS_intT myProc = globalClac.GetProcID();
  const FS_intT nProcs = globalClac.GetNProcs();

  std::vector<RawMessage> messages;
  messages.reserve(nProcs - 1);

  for(FS_intT sourceProc = 0; sourceProc < nProcs; ++sourceProc) {
    if(sourceProc == myProc)
      continue;

    globalClac.ReceiveBuffer(sourceProc);
    globalClac.SelectReceiveBuffer(sourceProc);

    RawMessage msg;
    globalClac.Unpack(&msg.worldProcID, 1);
    globalClac.Unpack(&msg.localProcID, 1);
    globalClac.Unpack(&msg.meshID, 1);
    FS_intT nFaces = 0;
    globalClac.Unpack(&nFaces, 1);

    for(FS_intT f = 0; f < nFaces; ++f) {
      RawFace raw;
      raw.originGlobalProc = msg.worldProcID;
      raw.meshID = msg.meshID;
      FS_intT nVerts = 0;
      globalClac.Unpack(&nVerts, 1);
      globalClac.Unpack(&raw.faceIndex, 1);
      raw.coords = FSFloatArrayT(nVerts, FS_3D);
      for(FS_intT v = 0; v < nVerts; ++v) {
        FS_float64T coords[3];
        globalClac.Unpack(coords, 3);
        for(FS_intT i = 0; i < FS_3D; i++)
          raw.coords(v, i) = coords[i];
      }
      raw.conn.Unpack(globalClac);
      msg.faces.push_back(std::move(raw));
    }
    messages.push_back(std::move(msg));
  }

  // 2 - reserve exact capacity for each group so connectivity never reallocates.
  AllGatheredFaces result;
  FS_intT nA = 0, nB = 0;
  for(const auto& msg : messages) {
    if(msg.meshID == 1)
      nA += static_cast<FS_intT>(msg.faces.size());
    else
      nB += static_cast<FS_intT>(msg.faces.size());
  }
  result.meshA.faces.connectivity.reserve(nA);
  result.meshA.faces.faces.reserve(nA);
  result.meshB.faces.connectivity.reserve(nB);
  result.meshB.faces.faces.reserve(nB);

  // 3 - Dispatch all the faces between mesh A and mesh B
  FS_intT nextIndexA = 0, nextIndexB = 0;

  for(auto& msg : messages) {

    const bool isMeshA = (msg.meshID == 1);
    GatheredFaces& g = isMeshA ? result.meshA : result.meshB;
    FS_intT& nextIndex = isMeshA ? nextIndexA : nextIndexB;
    // mapping local ProcID -> global ProcID
    g.localToGlobal[msg.localProcID] = msg.worldProcID;

    for(auto& raw : msg.faces) {
      g.cellToGlobalProc[{raw.conn.mOwner.mCellProcID, raw.conn.mOwner.mCell}] = raw.originGlobalProc;
      g.faces.connectivity.push_back(raw.conn);
      g.faces.faces.emplace_back(FSFace(g.faces.connectivity.back()), nextIndex++, raw.coords);
    }
  }

  return result;
}

} // namespace FSFaceExchange

namespace FSMatchExchange {

void Send(FSClac& clac, FS_intT destProc, const std::vector<FSFaceMatch>& matches)
{
  FSClac::sizeT bufSize = clac.GetBufSizeInt32(1); // n_matches
  for(const auto& m : matches)
    bufSize += m.GetBufSize(clac);

  clac.SelectSendBuffer(destProc);
  clac.InitSendBuffer(bufSize);
  FS_intT nm = static_cast<FS_intT>(matches.size());
  clac.Pack(&nm, 1);

  for(auto& m : matches) {
    FSFaceMatch fm = m;
    fm.Pack(clac);
  }

  clac.SendBuffer(destProc);
}

void ScatterSend(FSClac& clac, const std::vector<FSFaceMatch>& matches, const FSFaceExchange::GatheredFaces& gathered)
{
  std::map<FS_intT, std::vector<FSFaceMatch> > forAOrB;
  for(const auto& [local, global] : gathered.localToGlobal)
    forAOrB[global] = {};

  for(const auto& m : matches)
    forAOrB.at(gathered.cellToGlobalProc.at({m.ownerProc1, m.elemOwner1})).push_back(m);

  for(auto& [globalProc, subset] : forAOrB)
    Send(clac, globalProc, subset);
}

std::vector<FSFaceMatch> ScatterReceive(FSClac& clac, FS_intT clipperProc)
{

  clac.ReceiveBuffer(clipperProc);
  clac.SelectReceiveBuffer(clipperProc);
  FS_intT numMatches = 0;
  clac.Unpack(&numMatches, 1);

  std::vector<FSFaceMatch> matches(numMatches);
  for(FS_intT m = 0; m < numMatches; m++) {
    FSFaceMatch match;
    match.Unpack(clac);
    matches[m] = match;
  }
  return matches;
}

} // namespace FSMatchExchange

_FS_END_NAMESPACE
