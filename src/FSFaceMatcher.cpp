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
  // clippedPoly2D : taille + n * 2 coords
  s += clac.GetBufSizeInt32(1);
  s += clac.GetBufSizeFloat64(static_cast<FSClac::intT>(clippedPoly2D.size()) * 2);
  // clippedPoly3D : taille + n * 3 coords
  s += clac.GetBufSizeInt32(1);
  s += clac.GetBufSizeFloat64(static_cast<FSClac::intT>(clippedPoly3D.size()) * 3);
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
}

bool FSFaceMatcher::ComputeMatch(const FSClippingFace& f1,
                                 const FSClippingFace& f2,
                                 FSFaceMatch& out) const
{

  // 1. Identical faces
  if(FSClippingUtil::AreFacesIdentical(f1, f2, tol_)) {
    out.type = FSFaceMatch::IDENTICAL;
    return true;
  }

  // 2. Clipped algorithm
  auto poly1 = f1.projected2D();
  auto poly2 = FSClippingUtil::ProjectFaceInPlane(f2, f1);

  // We check that we have the same orientation for the two poly
  FSClippingUtil::EnsureSameOrientation(poly2, FSClippingUtil::PolygonSignedArea(poly1));

  std::vector<FSVec2> clip;
  if(!FSClippingUtil::ClipConvexPolygon(poly1, poly2, clip, tol_))
    return false;

  const FS_floatT areaClip = std::abs(FSClippingUtil::PolygonSignedArea(clip));
  const FS_floatT area2 = std::abs(FSClippingUtil::PolygonSignedArea(poly2));


  FSClippingUtil::RemoveDuplicatePoints(clip, 100 * tol_);
  if(clip.size() < 3)
    return false;

  if(FSClippingUtil::ArePointsColinear2D(clip, 100 * tol_))
    return false;

  if(areaClip < tol_)
    return false;

  // We do not test anymore for inclusion before the clipping algorithm because the inclusion clause is the source
  // of a lot of bugs (depends on orientation, the sign of the normal etc...).
  // All inclusions can be treated by the clipping algorithm without additionnal cost.
  // Finnaly we just add inclusion's match if the area of the poly and clipped faces are identical.
  // This remark was written after spending lots of time debugging the inclusion case.
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
  FSIntArrayT outIndicesBVHTree; // indice in the BVH Tree of the clipped faces
                                 // corresponding to the index in the
                                 // clippedFaces_=std::vector<FSClippingFace>

  for(const auto& subject : subjectFaces_) {

    FS_intT n = bvhClipped_.FindBoxesIntersectingWithBox(
      subject.boundingBox().boxMinMax, outIndicesBVHTree);
    FS_intT faceIndexSubject = subject.faceIndex();
    // here we can have potential face index clipped which are on other proc. We
    // need to acess all of them on this proc for now we'll supose that we are
    // on sequential
    outIndicesBVHTree.Resize(n);


    for(FS_intT k = 0; k < n; ++k) {
      FS_intT faceIndexClippedBVHTree = outIndicesBVHTree(k);
      FSFaceMatch match;
      match.face1 = faceIndexSubject;
      if(subject.topo()._faceFSDM) { // this is why we need the topo in FSClipping
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

      if(ComputeMatch(subject, clippedFaces_[faceIndexClippedBVHTree],
                      match)) {
        match.clippedPoly3D =
          FSClippingUtil::ProjectPoly2DTo3D(subject, match.clippedPoly2D);
        outMatches.emplace_back(std::move(match));
      }
    }
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
  }
}

_FS_END_NAMESPACE
