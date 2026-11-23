#ifndef FSVERTEXSNAP_H
#define FSVERTEXSNAP_H

#include "FSClipping/FSClippingFace.h"
#include <FSVec3.h>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

_FS_BEGIN_NAMESPACE

struct FSFaceMatch; // FSFaceMatcher.h owns the snapper

/*
  Merges the match vertices closer than snapTol, with one representative per side so no original node moves:
  an original node of that side, else the group point nearest to that side's face outline.
*/
class FSVertexSnap
{
public:
  // Groups are built once, on the raw matches (both frames) and the vertices of both face sets
  FSVertexSnap(const std::vector<FSClippingFace>& faces1, const std::vector<FSClippingFace>& faces2,
               const std::vector<FSFaceMatch>& rawMatches, FS_floatT snapTol);

  // side 0: clippedPoly3D in the frame of faces1, side 1: in the frame of faces2; returns the dropped pieces
  FS_intT Apply(std::vector<FSFaceMatch>& matches, FS_intT side) const;

private:
  using ExactKey = std::array<std::uint64_t, 3>;
  struct ExactKeyHash {
    std::size_t operator()(const ExactKey& k) const noexcept;
  };

  std::vector<FSVec3> points_;
  std::unordered_map<ExactKey, FS_intT, ExactKeyHash> indexOf_;
  std::array<std::vector<bool>, 2> isOriginal_;
  std::array<std::vector<FS_floatT>, 2> outlineDist_; // distance to the outline of its face, +inf if not on that side
  std::vector<FS_intT> group_;                        // point -> group
  std::array<std::vector<FS_intT>, 2> rep_;           // group -> representative point of each side, -1 if none

  FS_intT AddPoint(const FSVec3& p);
  FS_intT Find(const FSVec3& p) const;
};

_FS_END_NAMESPACE
#endif
