#ifndef FSFACESEPARATOR_H
#define FSFACESEPARATOR_H

#include "FSClac.h"
#include "FSCommon.h"
#include "FSError.h"
#include "FSLog.h"
#include "FSHashableFace.h"
#include "FSMesh.h"
#include "FSMeshCheck.h"
#include "FSMeshData.h"
#include "FSMeshFaceExtractor.h"
#include "FSMeshOpParams.h"
#include "FSMeshOpParamsArray.h"
#include "FSMeshPrintInfo.h"

#include "FSBoundaryFace.h"
#include "FSCellAdress.h"
#include "FSFace.h"
#include "FSClipping/FSBoundaryFaceProvider.h"

_FS_BEGIN_NAMESPACE

namespace FSFaceSeparator {

void SeparateFaces(const FSMesh& fsmesh, const FSMeshFaceExtractor& faceExtractor,
                   std::vector<FSBoundaryFace>& bdryFaces, const FS_intT markerBoundaryClipped);
} // namespace FSFaceSeparator
_FS_END_NAMESPACE

#endif