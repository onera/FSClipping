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
  an original node of that side, else a vertex lying on every face and face edge the group touches.
*/
class FSVertexSnap
{
public:
  // Built once on the raw matches (both frames) and the vertices of both face sets. Points closer than tol are one
  // vertex; vertices are then merged greedily by distance into groups that never hold two nodes of one mesh, stay
  // below snapTol in diameter and keep every point on its face outline
  FSVertexSnap(const std::vector<FSClippingFace>& faces1, const std::vector<FSClippingFace>& faces2,
               const std::vector<FSFaceMatch>& rawMatches, FS_floatT tol, FS_floatT snapTol);

  // side 0: clippedPoly3D in the frame of faces1, side 1: in the frame of faces2; returns the dropped pieces
  FS_intT Apply(std::vector<FSFaceMatch>& matches, FS_intT side) const;

  // Merges refused because the group would be too wide or would lose the outline on one side
  FS_intT NRefused() const noexcept { return nRefused_; }

private:
  using ExactKey = std::array<std::uint64_t, 3>;
  struct ExactKeyHash {
    std::size_t operator()(const ExactKey& k) const noexcept;
  };

  // Points closer than tol, with what they touch on each side
  struct FSSnapVertex {
    FS_intT point = -1; // its original node if any, else its first point
    std::array<bool, 2> original{false, false};
    std::array<FS_floatT, 2> outlineDist;           // distance to its face outline, +inf if not on that side
    std::array<std::vector<FS_intT>, 2> faces;      // sorted faces it belongs to
    std::array<std::vector<FS_intT>, 2> edges;      // sorted face edges it lies on
  };
  static constexpr FS_intT kNoRep = -1;     // the group has no valid representative on that side
  static constexpr FS_intT kNotOnSide = -2; // the group has no vertex on that side

  std::vector<FSVec3> points_;
  std::unordered_map<ExactKey, FS_intT, ExactKeyHash> indexOf_;
  std::vector<FSSnapVertex> vertices_;
  std::vector<FS_intT> group_;              // point -> group
  std::array<std::vector<FS_intT>, 2> rep_; // group -> representative point of each side, -1 if none
  FS_intT nRefused_ = 0;

  FS_intT Find(const FSVec3& p) const;
  // Representative vertex of a candidate group on side s, kNoRep or kNotOnSide
  FS_intT Representative(FS_intT s, const std::vector<FS_intT>& vertices) const;
};

_FS_END_NAMESPACE
#endif
