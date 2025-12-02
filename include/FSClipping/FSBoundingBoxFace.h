#ifndef FSBOUNDINGBOXFACE_H
#define FSBOUNDINGBOXFACE_H

#include "FSClac.h"
#include "FSCommon.h"
#include "FSConfig.h"
#include "FSFace.h"

_FS_BEGIN_NAMESPACE

struct FSBoundingBoxFace {
  // FSVec3 min;
  // FSVec3 max;
  FS_floatT boxMinMax[6];

  FSBoundingBoxFace()
      : boxMinMax(+FS_FLOATT_MAX, +FS_FLOATT_MAX, +FS_FLOATT_MAX,
                  -FS_FLOATT_MAX, -FS_FLOATT_MAX, -FS_FLOATT_MAX) {}

  void ExpandToInclude(const FSVec3 &p) {
    boxMinMax[0] = FSMin(boxMinMax[0], p[0]);
    boxMinMax[1] = FSMin(boxMinMax[1], p[1]);
    boxMinMax[2] = FSMin(boxMinMax[2], p[2]);
    boxMinMax[3] = FSMax(boxMinMax[3], p[0]);
    boxMinMax[4] = FSMax(boxMinMax[4], p[1]);
    boxMinMax[5] = FSMax(boxMinMax[5], p[2]);
  }
};

_FS_END_NAMESPACE

#endif