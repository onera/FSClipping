#include "FSError.h"

#include "FSClipping/FSClippingUtil.h"
#include "FSClipping/FSFaceMatcher.h"

_FS_BEGIN_NAMESPACE

bool FSFaceMatcher::ComputeMatch(const FSClippingFace& f1,
                                 const FSClippingFace& f2,
                                 FSFaceMatch& out) const
{

  // 1. identical test
  if(FSClippingUtil::AreFacesIdentical(f1, f2)) {
    out.type = FSFaceMatch::IDENTICAL;
    // out.intersectedArea = f1.projectedArea; // ou f2
    return true;
  }

  // 2. f2 projection in f1 plane
  std::vector<FSVec2> proj2 = FSClippingUtil::ProjectFaceInPlane(f2, f1);

  // 3. face inside test
  if(FSClippingUtil::IsFaceInsideTheOther(f1.projected2D(), proj2)) {
    out.type = FSFaceMatch::INCLUDED;
    // out.intersectedArea = PolygonArea(proj2);
    out.clippedPoly2D = proj2;
    return true;
  }

  // 4. inverse test
  std::vector<FSVec2> proj1 = FSClippingUtil::ProjectFaceInPlane(f1, f2);
  if(FSClippingUtil::IsFaceInsideTheOther(f2.projected2D(), proj1)) {
    out.type = FSFaceMatch::INCLUDED;
    // out.intersectedArea = PolygonArea(proj2);
    out.clippedPoly2D = proj1;
    return true;
  }

  // 4. Sutherland–Hodgman intersection
  std::vector<FSVec2> clip;
  if(!FSClippingUtil::ClipConvexPolygon(f1.projected2D(), proj2, clip))
    return false;

  out.type = FSFaceMatch::INTERSECTING;
  out.clippedPoly2D = std::move(clip);
  // out.intersectedArea = PolygonArea(clip);

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
    //   FSIntArrayT faceIndexClipped =
    //       bvhClipped_.GetKeysForCells(outIndicesBVHTree);

    for(FS_intT k = 0; k < n; ++k) {
      FS_intT faceIndexClippedBVHTree = outIndicesBVHTree(k);
      FSFaceMatch match;
      match.face1 = faceIndexSubject;
      match.elemOwner1 = subject.topo().GetOwnerCellFSDMIndex(); // this is why we need the topo in FSClipping
      match.elemOwnerType1 = subject.topo()._faceFSDM->mOwner.mCellType;
      // match.face2 = faceIndexClipped[k];
      match.face2 = clippedFaces_[faceIndexClippedBVHTree].faceIndex();
      match.elemOwner2 = clippedFaces_[faceIndexClippedBVHTree].topo().GetOwnerCellFSDMIndex();
      match.elemOwnerType1 = clippedFaces_[faceIndexClippedBVHTree].topo()._faceFSDM->mOwner.mCellType;
      // std::cout << faceIndexClipped[k] << "=?" << match.face2 << std::endl;
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
