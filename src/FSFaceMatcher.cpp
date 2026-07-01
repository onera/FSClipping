#include "FSClipping/FSClippingUtil.h"
#include "FSClipping/FSFaceMatcher.h"

_FS_BEGIN_NAMESPACE
FSClac::sizeT FSFaceMatch::GetBufSize(FSClac& clac) const
{
  FSClac::sizeT s = 0;
  // face1, face2, elemOwner1/2, faceOwner1/2, 4 CellTypes, type → 11 FS_intT
  s += clac.GetBufSizeInt32(11);
  // intersectedArea → 1 FS_float64T
  s += clac.GetBufSizeFloat64(1);
  // clippedPoly2D: size + n * 2 coords
  s += clac.GetBufSizeInt32(1);
  s += clac.GetBufSizeFloat64(static_cast<FSClac::intT>(clippedPoly2D.size()) * 2);
  // clippedPoly3D (face1 frame): size + n * 3 coords
  s += clac.GetBufSizeInt32(1);
  s += clac.GetBufSizeFloat64(static_cast<FSClac::intT>(clippedPoly3D.size()) * 3);
  // clippedPoly3D_face2 (face2 frame): size + n * 3 coords
  s += clac.GetBufSizeInt32(1);
  s += clac.GetBufSizeFloat64(static_cast<FSClac::intT>(clippedPoly3D_face2.size()) * 3);
  return s;
}

void FSFaceMatch::Pack(FSClac& clac)
{
  clac.Pack(&face1);
  clac.Pack(&face2);
  clac.Pack(&elemOwner1);
  clac.Pack(&elemOwner2);
  clac.Pack(&faceOwner1);
  clac.Pack(&faceOwner2);

  FS_intT t1 = elemOwnerType1, t2 = elemOwnerType2;
  FS_intT ft1 = faceOwnerType1, ft2 = faceOwnerType2;
  clac.Pack(&t1);
  clac.Pack(&t2);
  clac.Pack(&ft1);
  clac.Pack(&ft2);

  FS_intT matchType = static_cast<FS_intT>(type);
  clac.Pack(&matchType);

  clac.Pack(&intersectedArea);

  FS_intT n2 = static_cast<FS_intT>(clippedPoly2D.size());
  clac.Pack(&n2);
  for(auto& p : clippedPoly2D) {
    FS_float64T x = p[0], y = p[1];
    clac.Pack(&x);
    clac.Pack(&y);
  }

  FS_intT n3 = static_cast<FS_intT>(clippedPoly3D.size());
  clac.Pack(&n3);
  for(auto& p : clippedPoly3D) {
    FS_float64T x = p[0], y = p[1], z = p[2];
    clac.Pack(&x);
    clac.Pack(&y);
    clac.Pack(&z);
  }

  FS_intT n3f2 = static_cast<FS_intT>(clippedPoly3D_face2.size());
  clac.Pack(&n3f2);
  for(auto& p : clippedPoly3D_face2) {
    FS_float64T x = p[0], y = p[1], z = p[2];
    clac.Pack(&x);
    clac.Pack(&y);
    clac.Pack(&z);
  }
}

void FSFaceMatch::Unpack(FSClac& clac)
{
  clac.Unpack(&face1);
  clac.Unpack(&face2);
  clac.Unpack(&elemOwner1);
  clac.Unpack(&elemOwner2);
  clac.Unpack(&faceOwner1);
  clac.Unpack(&faceOwner2);

  FS_intT t1, t2, ft1, ft2;
  clac.Unpack(&t1);
  clac.Unpack(&t2);
  clac.Unpack(&ft1);
  clac.Unpack(&ft2);
  elemOwnerType1 = static_cast<FSMeshEnums::CellType>(t1);
  elemOwnerType2 = static_cast<FSMeshEnums::CellType>(t2);
  faceOwnerType1 = static_cast<FSMeshEnums::CellType>(ft1);
  faceOwnerType2 = static_cast<FSMeshEnums::CellType>(ft2);

  FS_intT matchType;
  clac.Unpack(&matchType);
  type = static_cast<MatchType>(matchType);

  clac.Unpack(&intersectedArea);

  FS_intT n2;
  clac.Unpack(&n2);
  clippedPoly2D.resize(n2);
  for(auto& p : clippedPoly2D) {
    FS_float64T x, y;
    clac.Unpack(&x);
    clac.Unpack(&y);
    p = FSVec2(x, y);
  }

  FS_intT n3;
  clac.Unpack(&n3);
  clippedPoly3D.resize(n3);
  for(auto& p : clippedPoly3D) {
    FS_float64T x, y, z;
    clac.Unpack(&x);
    clac.Unpack(&y);
    clac.Unpack(&z);
    p = FSVec3(x, y, z);
  }

  FS_intT n3f2;
  clac.Unpack(&n3f2);
  clippedPoly3D_face2.resize(n3f2);
  for(auto& p : clippedPoly3D_face2) {
    FS_float64T x, y, z;
    clac.Unpack(&x);
    clac.Unpack(&y);
    clac.Unpack(&z);
    p = FSVec3(x, y, z);
  }
}

bool FSFaceMatcher::ComputeMatch(const FSClippingFace& f1, const FSClippingFace& f2, FSFaceMatch& out) const
{
  // 1. Identical faces
  if(FSClippingUtil::AreFacesIdentical(f1, f2, tol_)) {
    out.type = FSFaceMatch::IDENTICAL;
    out.clippedPoly2D = f1.projected2D();
    return true;
  }

  // 2. Sutherland-Hodgman clipping in the plane of f1
  auto poly1 = f1.projected2D();
  auto poly2 = FSClippingUtil::ProjectFaceInPlane(f2, f1);

  FSClippingUtil::EnsureSameOrientation(poly2, FSClippingUtil::PolygonSignedArea(poly1));

  std::vector<FSVec2> clip;
  if(!FSClippingUtil::ClipConvexPolygon(poly1, poly2, clip, tol_))
    return false;

  const FS_floatT areaClip = std::abs(FSClippingUtil::PolygonSignedArea(clip));
  const FS_floatT area2 = std::abs(FSClippingUtil::PolygonSignedArea(poly2));

  FSClippingUtil::RemoveDuplicatePoints(clip, tol_);
  if(clip.size() < 3)
    return false;

  if(FSClippingUtil::ArePointsColinear2D(clip, tol_))
    return false;

  if(areaClip < tol_)
    return false;

  // Inclusion is detected via area equality after clipping — no dedicated
  // pre-check is needed. A separate inclusion branch was removed because it
  // was orientation-sensitive and produced incorrect results in edge cases.
  if(std::abs(areaClip - area2) < tol_)
    out.type = FSFaceMatch::INCLUDED;
  else
    out.type = FSFaceMatch::INTERSECTING;
  out.intersectedArea = areaClip;
  out.clippedPoly2D = std::move(clip);

  return true;
}

void FSFaceMatcher::ComputeMatches(std::vector<FSFaceMatch>& outMatches)
{
  outMatches.clear();
  FSIntArrayT outIndicesBVHTree; // indices in the BVH tree, parallel to clippedFaces_

  for(const auto& subject : subjectFaces_) {
    FS_intT n = bvhClipped_.FindBoxesIntersectingWithBox(subject.boundingBox().boxMinMax, outIndicesBVHTree);
    FS_intT faceIndexSubject = subject.faceIndex();
    outIndicesBVHTree.Resize(n);
    // FS_floatT areaAllClip = 0.0;
    // FS_floatT areaSubject = std::abs(FSClippingUtil::PolygonSignedArea(subject.projected2D()));

    for(FS_intT k = 0; k < n; ++k) {
      FS_intT faceIndexClippedBVHTree = outIndicesBVHTree(k);
      FSFaceMatch match;
      match.face1 = faceIndexSubject;
      if(subject.topo()._faceFSDM) {
        match.elemOwner1 = subject.topo().GetOwnerCellFSDMIndex();
        match.faceOwner1 = subject.topo().GetNeighborCellFSDMIndex();
        match.elemOwnerType1 = subject.topo()._faceFSDM->mOwner.mCellType;
        match.faceOwnerType1 = subject.topo()._faceFSDM->mNeighbor.mCellType;
      }

      match.face2 = clippedFaces_[faceIndexClippedBVHTree].faceIndex();
      if(clippedFaces_[faceIndexClippedBVHTree].topo()._faceFSDM) {
        match.elemOwner2 = clippedFaces_[faceIndexClippedBVHTree].topo().GetOwnerCellFSDMIndex();
        match.faceOwner2 = clippedFaces_[faceIndexClippedBVHTree].topo().GetNeighborCellFSDMIndex();
        match.elemOwnerType2 = clippedFaces_[faceIndexClippedBVHTree].topo()._faceFSDM->mOwner.mCellType;
        match.faceOwnerType2 = clippedFaces_[faceIndexClippedBVHTree].topo()._faceFSDM->mNeighbor.mCellType;
      }

      if(ComputeMatch(subject, clippedFaces_[faceIndexClippedBVHTree], match)) {
        match.clippedPoly3D = FSClippingUtil::ProjectPoly2DTo3D(subject, match.clippedPoly2D);
        // std::cout << "matche compute" << std::endl;
        //  Compute the intersection polygon in the clipped face's own frame by
        //  running the clipping algorithm with roles swapped. This ensures that
        //  vertices of the clipped face appear with their exact projected2D()
        //  values, making them consistent across all matches sharing that face.
        //  InvertMatches swaps clippedPoly3D <-> clippedPoly3D_face2 so that
        //  proc 1 always receives coordinates in its own face frame.
        FSFaceMatch match_swapped;
        if(ComputeMatch(clippedFaces_[faceIndexClippedBVHTree], subject, match_swapped)) {
          match.clippedPoly3D_face2 =
            FSClippingUtil::ProjectPoly2DTo3D(clippedFaces_[faceIndexClippedBVHTree], match_swapped.clippedPoly2D);
        }
        // areaAllClip += match.intersectedArea;
        outMatches.emplace_back(std::move(match));
      }
    }
    // if(std::abs(areaSubject - areaAllClip) < tol_) {
    //   std::cout << "You have a different area between the original face " << subject.faceIndex()
    //             << " and all the clipped " << std::endl;
    // }
  }
}

void FSFaceMatcher::InvertMatches(std::vector<FSFaceMatch>& matches)
{
  for(auto& m : matches) {
    std::swap(m.face1, m.face2);
    std::swap(m.elemOwner1, m.elemOwner2);
    std::swap(m.elemOwnerType1, m.elemOwnerType2);
    std::swap(m.faceOwner1, m.faceOwner2);
    std::swap(m.faceOwnerType1, m.faceOwnerType2);
    std::swap(m.clippedPoly3D, m.clippedPoly3D_face2);
  }
}

_FS_END_NAMESPACE
