#include "FSClipping/FSFaceExchange.h"
#include <unordered_map>

_FS_BEGIN_NAMESPACE

namespace FSFaceExchange {

void Send(FSClac& clac, FS_intT destProc, const std::vector<FSClippingFace>& faces)
{
  // Per face: [n_verts (int32)] [coords (float64 * n*3)] [owner (5 int32)] [neighbor (5 int32)]
  FSClac::sizeT bufSize = clac.GetBufSizeInt32(1); // n_faces
  for(const auto& f : faces) {
    bufSize += clac.GetBufSizeInt32(1);                                                    // n_verts
    bufSize += clac.GetBufSizeInt32(1);                                                    // faceIndex
    bufSize += clac.GetBufSizeFloat64(static_cast<FSClac::intT>(f.vertices().size()) * 3); // coords
    bufSize += f.topo()._faceFSDM->GetBufSize(clac);
  }

  clac.SelectSendBuffer(destProc);
  clac.InitSendBuffer(bufSize);

  FS_intT n = static_cast<FS_intT>(faces.size());
  clac.Pack(&n, 1);

  for(const auto& f : faces) {
    FS_intT nv = static_cast<FS_intT>(f.vertices().size());
    clac.Pack(&nv, 1);
    FS_intT faceIndex = f.faceIndex();
    clac.Pack(&faceIndex, 1);
    for(const auto& v : f.vertices()) {
      FS_float64T coords[3] = {v[0], v[1], v[2]};
      clac.Pack(coords, 3);
    }
    // Pack FSDM connectivity (owner + neighbor cell face, 5 FS_intT each)
    FSCellFace owner = f.topo()._faceFSDM->mOwner;
    FSCellFace neighbor = f.topo()._faceFSDM->mNeighbor;
    owner.Pack(clac);
    neighbor.Pack(clac);
  }

  clac.SendBuffer(destProc);
}

void GatherSend(FSClac& clac, FS_intT clipperProc, FS_intT meshID, const std::vector<FSClippingFace>& faces)
{
  FS_intT worldProcID = clac.GetWorldProcID();
  FS_intT localProcID = clac.GetProcID();

  // Header: worldProcID + localProcID + meshID + n_faces
  FSClac::sizeT bufSize = clac.GetBufSizeInt32(4);
  for(const auto& f : faces) {
    bufSize += clac.GetBufSizeInt32(2);
    bufSize += clac.GetBufSizeFloat64(static_cast<FSClac::intT>(f.vertices().size()) * 3);
    bufSize += f.topo()._faceFSDM->GetBufSize(clac);
  }

  clac.SelectSendBuffer(clipperProc);
  clac.InitSendBuffer(bufSize);
  clac.Pack(&worldProcID, 1);
  clac.Pack(&localProcID, 1);
  clac.Pack(&meshID, 1);
  FS_intT n = static_cast<FS_intT>(faces.size());
  clac.Pack(&n, 1);

  for(const auto& f : faces) {
    FS_intT nv = static_cast<FS_intT>(f.vertices().size());
    clac.Pack(&nv, 1);
    FS_intT faceIndex = f.faceIndex();
    clac.Pack(&faceIndex, 1);
    for(const auto& v : f.vertices()) {
      FS_float64T coords[3] = {v[0], v[1], v[2]};
      clac.Pack(coords, 3);
    }
    FSCellFace owner = f.topo()._faceFSDM->mOwner;
    FSCellFace neighbor = f.topo()._faceFSDM->mNeighbor;
    owner.Pack(clac);
    neighbor.Pack(clac);
  }

  clac.SendBuffer(clipperProc);
}

AllGatheredFaces GatherReceiveAll(FSClac& clac)
{
  // Pass 1 — receive raw data from every non-clipper proc.
  // Raw storage prevents connectivity vector reallocation (which would dangle FSFace pointers).
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

  const FS_intT myProc = clac.GetProcID();
  const FS_intT nProcs = clac.GetNProcs();

  std::vector<RawMessage> messages;
  messages.reserve(nProcs - 1);

  for(FS_intT sourceProc = 0; sourceProc < nProcs; ++sourceProc) {
    if(sourceProc == myProc)
      continue;

    clac.ReceiveBuffer(sourceProc);
    clac.SelectReceiveBuffer(sourceProc);

    RawMessage msg;
    clac.Unpack(&msg.worldProcID, 1);
    clac.Unpack(&msg.localProcID, 1);
    clac.Unpack(&msg.meshID, 1);
    FS_intT nFaces = 0;
    clac.Unpack(&nFaces, 1);

    for(FS_intT f = 0; f < nFaces; ++f) {
      RawFace raw;
      raw.originGlobalProc = msg.worldProcID;
      raw.meshID = msg.meshID;
      FS_intT nVerts = 0;
      clac.Unpack(&nVerts, 1);
      clac.Unpack(&raw.faceIndex, 1);
      raw.coords = FSFloatArrayT(nVerts, FS_3D);
      for(FS_intT v = 0; v < nVerts; ++v) {
        FS_float64T coords[3];
        clac.Unpack(coords, 3);
        for(FS_intT i = 0; i < FS_3D; i++)
          raw.coords(v, i) = coords[i];
      }
      raw.conn.Unpack(clac);
      msg.faces.push_back(std::move(raw));
    }
    messages.push_back(std::move(msg));
  }

  // Pass 2 — reserve exact capacity for each group so connectivity never reallocates.
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

  // Face indices are local to the sending proc, so faces coming from different
  // procs of one mesh collide on the same index. Downstream those indices are the
  // identity of a face (FSFaceMatch::face1/face2, and the per-face coverage sum in
  // PreserveUncoveredFaces), so a collision made two distinct faces look like one:
  // their intersection areas were summed, the total came to twice the face area,
  // and the face was wrongly reported as not covered. Renumber the concatenated
  // list contiguously so each gathered face has a unique index on the clipper.
  FS_intT nextIndexA = 0, nextIndexB = 0;

  for(auto& msg : messages) {
    const bool isMeshA = (msg.meshID == 1);
    GatheredFaces& g = isMeshA ? result.meshA : result.meshB;
    FS_intT& nextIndex = isMeshA ? nextIndexA : nextIndexB;
    g.localToGlobal[msg.localProcID] = msg.worldProcID;
    for(auto& raw : msg.faces) {
      // Key on (owner proc, cell): cell IDs are local, so two procs of this same
      // mesh can both own a cell numbered e.g. 201.
      g.cellToGlobalProc[{raw.conn.mOwner.mCellProcID, raw.conn.mOwner.mCell}] = raw.originGlobalProc;
      g.faces.connectivity.push_back(raw.conn);
      g.faces.faces.emplace_back(FSFace(g.faces.connectivity.back()), nextIndex++, raw.coords);
    }
  }

  return result;
}

ReceivedFaces Receive(FSClac& clac, FS_intT sourceProc)
{
  clac.ReceiveBuffer(sourceProc);
  clac.SelectReceiveBuffer(sourceProc);

  ReceivedFaces result;

  FS_intT nFaces = 0;
  clac.Unpack(&nFaces, 1);
  result.faces.reserve(nFaces);
  result.connectivity.reserve(nFaces);

  for(FS_intT f = 0; f < nFaces; ++f) {
    FS_intT nVerts = 0;
    clac.Unpack(&nVerts, 1);
    FS_intT faceIndex = 0;
    clac.Unpack(&faceIndex, 1);

    FSFloatArrayT faceNodeCoordinates(nVerts, FS_3D);
    for(FS_intT v = 0; v < nVerts; ++v) {
      FS_float64T coords[3];
      clac.Unpack(coords, 3);
      for(FS_intT i = 0; i < FS_3D; i++)
        faceNodeCoordinates(v, i) = coords[i];
    }

    FSFaceConnectivity conn;
    conn.Unpack(clac);

    result.connectivity.push_back(conn);
    result.faces.emplace_back(FSFace(result.connectivity.back()), faceIndex, faceNodeCoordinates);
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

std::vector<FSFaceMatch> Receive(FSClac& clac, FS_intT sourceProc)
{
  clac.ReceiveBuffer(sourceProc);
  clac.SelectReceiveBuffer(sourceProc);
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

static FSFaceMatch InvertSingle(const FSFaceMatch& m)
{
  FSFaceMatch inv = m;
  std::swap(inv.face1, inv.face2);
  std::swap(inv.elemOwner1, inv.elemOwner2);
  std::swap(inv.ownerProc1, inv.ownerProc2);
  std::swap(inv.elemOwnerType1, inv.elemOwnerType2);
  std::swap(inv.faceOwner1, inv.faceOwner2);
  std::swap(inv.faceOwnerType1, inv.faceOwnerType2);
  std::swap(inv.clippedPoly3D, inv.clippedPoly3D_face2);
  return inv;
}

void ScatterSend(FSClac& clac, const std::vector<FSFaceMatch>& matches, const FSFaceExchange::GatheredFaces& gathered)
{
  // Initialize empty buckets for every mesh proc so procs with 0 matches
  // still receive a message and ScatterReceive never hangs.
  std::map<FS_intT, std::vector<FSFaceMatch> > forAOrB;
  for(const auto& [local, global] : gathered.localToGlobal)
    forAOrB[global] = {};

  // Route on (owner proc, cell), the only globally unique cell identity: the
  // several procs of this mesh may each number a different cell with the same ID.
  for(const auto& m : matches)
    forAOrB.at(gathered.cellToGlobalProc.at({m.ownerProc1, m.elemOwner1})).push_back(m);

  for(auto& [globalProc, subset] : forAOrB)
    Send(clac, globalProc, subset);
}

std::vector<FSFaceMatch> ScatterReceive(FSClac& clac, FS_intT clipperProc) { return Receive(clac, clipperProc); }

} // namespace FSMatchExchange

_FS_END_NAMESPACE
