#include "FSClipping/FSFaceExchange.h"

_FS_BEGIN_NAMESPACE

namespace FSFaceExchange {

void Send(FSClac& clac, FS_intT destProc, const std::vector<FSClippingFace>& faces)
{
  FSClac::sizeT bufSize = clac.GetBufSizeInt32(1); // n_faces
  for(const auto& f : faces)
    bufSize += clac.GetBufSizeInt32(1) + clac.GetBufSizeFloat64(static_cast<FSClac::intT>(f.vertices().size()) * 3);

  clac.SelectSendBuffer(destProc);
  clac.InitSendBuffer(bufSize);

  FS_intT n = static_cast<FS_intT>(faces.size());
  clac.Pack(&n, 1);
  for(const auto& f : faces) {
    FS_intT nv = static_cast<FS_intT>(f.vertices().size());
    clac.Pack(&nv, 1);
    for(const auto& v : f.vertices()) {
      FS_float64T coords[3] = {v[0], v[1], v[2]};
      clac.Pack(coords, 3);
    }
  }

  clac.SendBuffer(destProc);
}

std::vector<FSClippingFace> Receive(FSClac& clac, FS_intT sourceProc)
{
  clac.ReceiveBuffer(sourceProc);
  clac.SelectReceiveBuffer(sourceProc);

  FS_intT nFaces = 0;
  clac.Unpack(&nFaces, 1);
  std::vector<FSClippingFace> faces;
  faces.reserve(nFaces);
  for(FS_intT f = 0; f < nFaces; ++f) {
    FS_intT nVerts = 0;
    clac.Unpack(&nVerts, 1);
    std::vector<FSVec3> verts(nVerts);
    for(FS_intT v = 0; v < nVerts; ++v) {
      FS_float64T coords[3];
      clac.Unpack(coords, 3);
      verts[v] = FSVec3(coords[0], coords[1], coords[2]);
    }
    faces.emplace_back(verts);
  }
  return faces;
}

} // namespace FSFaceExchange

_FS_END_NAMESPACE
