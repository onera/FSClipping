#include "FSClipping/FSFaceExchange.h"
#include <unordered_map>

_FS_BEGIN_NAMESPACE

namespace FSFaceExchange {

void Send(FSClac& clac, FS_intT destProc, const std::vector<FSClippingFace>& faces)
{
  // --- Compute exact buffer size ---
  // Per face: [n_verts (int32)] [coords (float64 * n*3)] [owner (5 int32)] [neighbor (5 int32)]
  FSClac::sizeT bufSize = clac.GetBufSizeInt32(1); // n_faces
  for(const auto& f : faces) {
    bufSize += clac.GetBufSizeInt32(1);                                                            // n_verts
    bufSize += clac.GetBufSizeFloat64(static_cast<FSClac::intT>(f.vertices().size()) * 3);        // coords
    bufSize += clac.GetBufSizeInt32(10); // FSCellFace owner (5 int32) + neighbor (5 int32)
  }

  clac.SelectSendBuffer(destProc);
  clac.InitSendBuffer(bufSize);

  // --- Pack ---
  FS_intT n = static_cast<FS_intT>(faces.size());
  clac.Pack(&n, 1);

  for(const auto& f : faces) {
    // geometry
    FS_intT nv = static_cast<FS_intT>(f.vertices().size());
    clac.Pack(&nv, 1);
    for(const auto& v : f.vertices()) {
      FS_float64T coords[3] = {v[0], v[1], v[2]};
      clac.Pack(coords, 3);
    }
    // FSDM connectivity (owner + neighbor FSCellFace, 5 FS_intT each)
    FSCellFace owner   = f.topo()._faceFSDM->mOwner;
    FSCellFace neighbor = f.topo()._faceFSDM->mNeighbor;
    owner.Pack(clac);
    neighbor.Pack(clac);
  }

  clac.SendBuffer(destProc);
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
    // geometry
    FS_intT nVerts = 0;
    clac.Unpack(&nVerts, 1);
    std::vector<FSVec3> verts(nVerts);
    for(FS_intT v = 0; v < nVerts; ++v) {
      FS_float64T coords[3];
      clac.Unpack(coords, 3);
      verts[v] = FSVec3(coords[0], coords[1], coords[2]);
    }
    result.faces.emplace_back(verts);

    // FSDM connectivity
    FSFaceConnectivity conn;
    conn.mOwner.Unpack(clac);
    conn.mNeighbor.Unpack(clac);
    result.connectivity.push_back(conn);
  }

  return result;
}

void InjectTopology(std::vector<FSFaceMatch>& matches, const ReceivedFaces& subjectReceived)
{
  // Build a map: faceIndex() → position in received.connectivity
  std::unordered_map<FS_intT, FS_intT> faceIdxToConn;
  faceIdxToConn.reserve(subjectReceived.faces.size());
  for(FS_intT i = 0; i < static_cast<FS_intT>(subjectReceived.faces.size()); ++i)
    faceIdxToConn[subjectReceived.faces[i].faceIndex()] = i;

  for(auto& match : matches) {
    auto it = faceIdxToConn.find(match.face1);
    if(it == faceIdxToConn.end())
      continue;
    const auto& conn = subjectReceived.connectivity[it->second];
    match.elemOwner1    = conn.mOwner.mCell;
    match.elemOwnerType1 = conn.mOwner.mCellType;
    match.faceOwner1    = conn.mNeighbor.mCell;
    match.faceOwnerType1 = conn.mNeighbor.mCellType;
  }
}

} // namespace FSFaceExchange

_FS_END_NAMESPACE
