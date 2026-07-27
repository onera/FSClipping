#include "FSClipping/FSClippingUtil.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSCell2NodeBuilder.h" // NodeKey

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

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
  // nodeGlobalIds: size + n ints ; globalCellId: 1 int
  s += clac.GetBufSizeInt32(1);
  s += clac.GetBufSizeInt32(static_cast<FSClac::intT>(nodeGlobalIds.size()));
  s += clac.GetBufSizeInt32(1);
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

  FS_intT nIds = static_cast<FS_intT>(nodeGlobalIds.size());
  clac.Pack(&nIds);
  for(auto& id : nodeGlobalIds)
    clac.Pack(&id);
  clac.Pack(&globalCellId);
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

  FS_intT nIds;
  clac.Unpack(&nIds);
  nodeGlobalIds.resize(nIds);
  for(auto& id : nodeGlobalIds)
    clac.Unpack(&id);
  clac.Unpack(&globalCellId);
}

bool FSFaceMatcher::ComputeMatch(const FSClippingFace& f1, const FSClippingFace& f2, FSFaceMatch& out) const
{
  // 1. Identical faces
  if(FSClippingUtil::AreFacesIdentical(f1, f2, tol_)) {
    out.type = FSFaceMatch::IDENTICAL;
    out.clippedPoly2D = f1.projected2D();
    out.intersectedArea = std::abs(FSClippingUtil::PolygonSignedArea(out.clippedPoly2D));
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

  // Mean width of the clip polygon (2*area/perimeter) below tol means the
  // intersection is a sliver of negligible thickness. Comparing the raw area
  // (length^2) against tol (length) was scale-dependent and rejected real
  // intersections on small meshes.
  const FS_floatT perimClip = FSClippingUtil::PolygonPerimeter(clip);
  if(2.0 * areaClip < tol_ * perimClip)
    return false;

  // Inclusion is detected via area equality after clipping — no dedicated
  // pre-check is needed. A separate inclusion branch was removed because it
  // was orientation-sensitive and produced incorrect results in edge cases.
  // The missing strip between clip and poly2 has area ~= width * perimeter/2,
  // so the threshold is expressed with the same dimensions.
  if(std::abs(areaClip - area2) < tol_ * 0.5 * FSClippingUtil::PolygonPerimeter(poly2))
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

  // Expand query bbox by tol_ to account for floating-point imprecision at flat
  // interfaces (e.g. two meshes whose shared plane coordinate differs by ~1 ULP).
  FS_floatT expandedBox[6];
  for(const auto& subject : subjectFaces_) {
    const FS_floatT* rawBox = subject.boundingBox().boxMinMax;
    for(FS_intT d = 0; d < FS_3D; ++d) {
      expandedBox[d] = rawBox[d] - tol_;
      expandedBox[d + FS_3D] = rawBox[d + FS_3D] + tol_;
    }
    FS_intT n = bvhClipped_.FindBoxesIntersectingWithBox(expandedBox, outIndicesBVHTree);
    FS_intT faceIndexSubject = subject.faceIndex();
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
        // Compute the intersection polygon in the clipped face's own frame by
        // running the clipping algorithm with roles swapped. This ensures that
        // vertices of the clipped face appear with their exact projected2D()
        // values, making them consistent across all matches sharing that face.
        // InvertMatches swaps clippedPoly3D <-> clippedPoly3D_face2 so that
        // proc 1 always receives coordinates in its own face frame.
        FSFaceMatch match_swapped;
        if(ComputeMatch(clippedFaces_[faceIndexClippedBVHTree], subject, match_swapped)) {
          match.clippedPoly3D_face2 =
            FSClippingUtil::ProjectPoly2DTo3D(clippedFaces_[faceIndexClippedBVHTree], match_swapped.clippedPoly2D);
        }
        outMatches.emplace_back(std::move(match));
      }
    }
  }
  // Keep the raw intersections so that ComputeInvertedMatches can apply the
  // coverage criterion independently on the clipped side.
  rawMatches_ = outMatches;

  PreserveUncoveredFaces(outMatches, subjectFaces_);
}

void FSFaceMatcher::PreserveUncoveredFaces(std::vector<FSFaceMatch>& matches,
                                           const std::vector<FSClippingFace>& faces) const
{
  // Sum of the intersection areas per face (the face is identified by face1)
  std::unordered_map<FS_intT, FS_floatT> coveredArea;
  for(const auto& m : matches)
    coveredArea[m.face1] += m.intersectedArea;

  // A face whose clipped polygons do not cover its whole area would leave a
  // hole in the rebuilt surface (e.g. rim faces of two cylinders in relative
  // rotation, only partially covered by the other mesh).
  std::unordered_set<FS_intT> uncovered;
  for(const auto& face : faces) {
    const FS_floatT faceArea = std::abs(FSClippingUtil::PolygonSignedArea(face.projected2D()));
    const auto it = coveredArea.find(face.faceIndex());
    const FS_floatT covered = (it == coveredArea.end()) ? 0.0 : it->second;

    // An uncovered strip of width tol along the face boundary has area
    // ~= tol * perimeter/2 — dimensionally consistent threshold (see ComputeMatch).
    // const FS_floatT perimeter = FSClippingUtil::PolygonPerimeter(face.projected2D());
    if(std::abs(covered - faceArea) > tol_)
      uncovered.insert(face.faceIndex());
  }

  if(uncovered.empty())
    return;

  // Drop the partial clips of the uncovered faces...
  matches.erase(std::remove_if(matches.begin(), matches.end(),
                               [&uncovered](const FSFaceMatch& m) { return uncovered.contains(m.face1); }),
                matches.end());

  // ...and keep each uncovered face whole, as a single NOT_COVERED match
  // carrying the original face polygon.
  for(const auto& face : faces) {
    if(!uncovered.contains(face.faceIndex()))
      continue;

    FSFaceMatch match;
    match.face1 = face.faceIndex();
    if(face.topo()._faceFSDM) {
      match.elemOwner1 = face.topo().GetOwnerCellFSDMIndex();
      match.faceOwner1 = face.topo().GetNeighborCellFSDMIndex();
      match.elemOwnerType1 = face.topo()._faceFSDM->mOwner.mCellType;
      match.faceOwnerType1 = face.topo()._faceFSDM->mNeighbor.mCellType;
    }
    match.type = FSFaceMatch::NOT_COVERED;
    match.intersectedArea = std::abs(FSClippingUtil::PolygonSignedArea(face.projected2D()));
    match.clippedPoly2D = face.projected2D();
    match.clippedPoly3D = face.vertices();
    matches.emplace_back(std::move(match));
  }
}

FS_intT FSFaceMatcher::AssignGlobalNodeIds(std::vector<FSFaceMatch>& matches, FS_floatT tol)
{
  std::unordered_map<NodeKey, FS_intT> nodeIds;
  FS_intT nextId = 0;

  for(std::size_t m = 0; m < matches.size(); ++m) {
    auto& match = matches[m];
    match.globalCellId = static_cast<FS_intT>(m);

    match.nodeGlobalIds.resize(match.clippedPoly3D.size());
    for(std::size_t i = 0; i < match.clippedPoly3D.size(); ++i) {
      NodeKey key(match.clippedPoly3D[i], tol);
      auto [it, inserted] = nodeIds.try_emplace(key, nextId);
      if(inserted)
        ++nextId;
      match.nodeGlobalIds[i] = it->second;
    }
  }

  return nextId;
}

void FSFaceMatcher::ComputeInvertedMatches(std::vector<FSFaceMatch>& outMatches) const
{
  // Start from the raw matches, not from the subject-side list: a partial clip
  // dropped because its subject face is NOT_COVERED is still a valid
  // intersection for the clipped face it belongs to.
  outMatches = rawMatches_;

  for(auto& m : outMatches) {
    std::swap(m.face1, m.face2);
    std::swap(m.elemOwner1, m.elemOwner2);
    std::swap(m.elemOwnerType1, m.elemOwnerType2);
    std::swap(m.faceOwner1, m.faceOwner2);
    std::swap(m.faceOwnerType1, m.faceOwnerType2);
    std::swap(m.clippedPoly3D, m.clippedPoly3D_face2);
  }

  PreserveUncoveredFaces(outMatches, clippedFaces_);
}

_FS_END_NAMESPACE
