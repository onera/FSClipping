#ifndef FACEFSDM_HPP
#define FACEFSDM_HPP

#include "FSMesh/FSHashableFace.h"

/// Structure to store FSDM info on a one to one connectivity through a face.
struct FSFace {
  /// The associated face in the FSMeshFaceExtractor, always non-null (but not a
  /// reference, so we can copy / swap).
  const FSFaceConnectivity *_faceFSDM;

  /// Constructor.
  /**
   * @param face The associated face in the FSMeshFaceExtractor to point on.
   */
  explicit FSFace(const FSFaceConnectivity &face) : _faceFSDM(&face) {}

  /// Get the FSDM element link owner cell index.
  FS_intT GetOwnerCellFSDMIndex() const { return _faceFSDM->mOwner.mCell; }

  /// Get the FSDM element link neighbor cell index.
  FS_intT GetNeighborCellFSDMIndex() const {
    return _faceFSDM->mNeighbor.mCell;
  }

  /// Get the FSDM element link neighbor cell process index.
  FS_intT GetNeighborCellProcFSDMIndex() const {
    return _faceFSDM->mNeighbor.mCellProcID;
  }

  /// Is this less than other?
  /**
   * @param other The object to compare against.
   * @return true if this object is ordered before other, false otherwise.
   */
  bool operator<(const FSFace &other) const {
    return *_faceFSDM < *other._faceFSDM;
  }

  // --- Default constructor ---
  FSFace() : _faceFSDM() {
    FSFaceConnectivity faceConnectivityDummy;
    this->_faceFSDM = &faceConnectivityDummy;
  }
};

#endif