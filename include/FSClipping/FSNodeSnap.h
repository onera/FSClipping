#ifndef FSNODESNAP_H
#define FSNODESNAP_H

#include "FSClipping/FSClippingFace.h"
#include <FSVec3.h>
#include <vector>

_FS_BEGIN_NAMESPACE

// Displacement of one interface node, identified by its exact original coordinates
struct FSNodeMove {
  FSVec3 from;
  FSVec3 to;
};

// Snaps the interface nodes of one mesh onto the nearly coincident nodes of the other one
class FSNodeSnap
{
public:
  // Moves of the vertices of `moving` lying within snapTol of a vertex of `target`
  static std::vector<FSNodeMove> Compute(const std::vector<FSClippingFace>& target,
                                         const std::vector<FSClippingFace>& moving, FS_floatT snapTol);

  static void Apply(const std::vector<FSNodeMove>& moves, std::vector<FSClippingFace>& faces);

  // Returns the number of moved points
  static FS_intT Apply(const std::vector<FSNodeMove>& moves, std::vector<FSVec3>& points);
};

_FS_END_NAMESPACE
#endif
