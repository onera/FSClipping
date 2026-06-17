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

} // namespace FSMatchExchange

_FS_END_NAMESPACE
