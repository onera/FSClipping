#ifndef FSFACESEPARATOR_H
#define FSFACESEPARATOR_H

#include "FSClac.h"
#include "FSCommon.h"
#include "FSError.h"
#include "FSLog.h"
#include "FSMesh/FSHashableFace.h"
#include "FSMesh/FSMesh.h"
#include "FSMesh/FSMeshCheck.h"
#include "FSMesh/FSMeshData.h"
#include "FSMesh/FSMeshFaceExtractor.h"
#include "FSMesh/FSMeshOpParams.h"
#include "FSMesh/FSMeshOpParamsArray.h"
#include "FSMesh/FSMeshPrintInfo.h"

#include "FSBoundaryFace.h"
#include "FSCellAdress.h"
#include "FSFace.h"

_FS_BEGIN_NAMESPACE

namespace FSFaceSeparator {

void SeparateFaces(const FSMesh &fsmesh,
                   const FSMeshFaceExtractor &faceExtractor,
                   std::vector<FSFace> &localFacesFSDM,
                   std::vector<FSFace> &haloFacesFSDM,
                   std::vector<FSBoundaryFace> &bdryFacesFSDM,
                   std::vector<FSFace> &extBdryFacesFSDM,
                   std::vector<FSBoundaryFace> &sepBdryFacesFSDM);

bool PrepareNodalData(const FSMesh &fsmesh, FSMeshFaceExtractor &faceExtractor);

std::vector<FSBoundaryFace> SeparateBoundariesFaceWithMarker(
    FSMesh &fsmesh, FSMeshFaceExtractor &faceExtractor, const FS_intT marker);

} // namespace FSFaceSeparator
_FS_END_NAMESPACE

#endif