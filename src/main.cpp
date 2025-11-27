#include <cstdlib>
#include <iostream>

#include "FSClipping/main.h"
#include "FSCommon.h"
#include "FSVec3.h"
#include "FSConfig.h"

_FS_BEGIN_NAMESPACE

int main()
{
    FSVec3 u;
    std::cout << u.L1Norm() << std::endl;

    return EXIT_SUCCESS;
}

_FS_END_NAMESPACE