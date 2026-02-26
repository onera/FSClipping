#include "FSError.h"

#include "FSClipping/FSClippingUtil.h"
#include "FSClipping/FSFaceMatcher.h"

_FS_BEGIN_NAMESPACE

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

  // We do not test anymore for inclusion before the clipping algorithm because the inclusion clause is the source
  // of a lot of bugs (depends on orientation, the sign of the normal etc...).
  // All inclusions can be treated by the clipping algorithm without additionnal cost.
  // Finnaly we just add inclusion's match if the area of the poly and clipped faces are identical.
  // This remark was written after spending lots of time debugging the inclusion case.
  if(std::abs(areaClip - area2) < tol_)
    out.type = FSFaceMatch::INCLUDED;
  else
    out.type = FSFaceMatch::INTERSECTING;

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
      match.elemOwner1 = subject.topo().GetOwnerCellFSDMIndex(); // this is why we need the topo in FSClipping
      match.elemOwnerType1 = subject.topo()._faceFSDM->mOwner.mCellType;
      match.face2 = clippedFaces_[faceIndexClippedBVHTree].faceIndex();
      match.elemOwner2 = clippedFaces_[faceIndexClippedBVHTree].topo().GetOwnerCellFSDMIndex();
      match.elemOwnerType2 = clippedFaces_[faceIndexClippedBVHTree].topo()._faceFSDM->mOwner.mCellType;
      if(ComputeMatch(subject, clippedFaces_[faceIndexClippedBVHTree],
                      match)) {
        match.clippedPoly3D =
          FSClippingUtil::ProjectPoly2DTo3D(subject, match.clippedPoly2D);
        outMatches.emplace_back(std::move(match));
      }
    }
  }
}

_FS_END_NAMESPACE
