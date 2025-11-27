#include "FSCommon.h"
#include <FSVec3.h>

#include "gtest/gtest.h"


_FS_BEGIN_NAMESPACE

TEST(FSSandBox, FSSandBox)
{
    FSVec3 u(1.0, 2.3, 4.3);
    ASSERT_TRUE(u.Norm() > 0);
}

_FS_END_NAMESPACE