#ifndef BDRYFACEFSDM_H
#define BDRYFACEFSDM_H

#include "FSFace.h"

_FS_BEGIN_NAMESPACE

/// Structure to store info on a local FSDM boundary face.
struct FSBoundaryFace : FSFace {
  using BaseClass = FSFace;

  /// FSDM boundary marker, which this boundary face belongs to.
  FS_intT _marker;
  /// FSDM boundary face index uniquely defining this face in the global mesh.
  FS_intT _fsdmFaceGlobalNumber;

  /// Constructor.
  /**
   * @param face The associated face in the FSMeshFaceExtractor to point on.
   * @param marker The FSDM boundary marker, which this boundary face belongs
   * to.
   * @param fsdmFaceGlobalNumber The FSDM boundary face index uniquely defining
   * this face in the global mesh.
   */
  FSBoundaryFace(const FSFaceConnectivity& face, const FS_intT marker, const FS_intT fsdmFaceGlobalNumber)
    : FSFace(face), _marker(marker), _fsdmFaceGlobalNumber(fsdmFaceGlobalNumber)
  {
  }

  _FS_END_NAMESPACE
};

#endif