import FSDM
import os, math
from FSDataManager import FSClac, FSLog, FSError, FSMesh, FSMeshEnums, FSQuantityDescArray, FSUtil, FSFloatArray, FSStringArray
from FSDataManager import FSDataName, FSDataSpecArray, FSDatasetInfo, FSQuantityDesc, FSDataValue
from FSDataManager import FSDataManager

import FSClipping

# --- instantiate clac and data manager
clac = FSClac()
dm = FSDataManager(clac)

# --- instantiate and register FSMesh objects
fsmeshOrig1 = dm.GetMesh("original1")
fsmeshOrig2 = dm.GetMesh("original2")
fsmeshClipped = dm.GetMesh("clippedMesh", True)

meshOps = (("ImportMeshTAU", {"MeshFilename"    : "/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/input/cube_hexa_coarse.grid",}),
           "PrintInfo",
           "Check",
          )
fsmeshOrig1.DoOps(meshOps) or FSError.PrintAndExit()

meshOps = (("ImportMeshTAU", {"MeshFilename"    : "/stck/aleprevo/code_dev/CODA_src/FSClipping/test/Mesh/input/cube_tetra_fine.grid",}),
           "PrintInfo",
           "Check",
          )

fsmeshOrig2.DoOps(meshOps) or FSError.PrintAndExit()

dataManagerOps = (("ClippedMesh", {"MeshKeyOrig1"   : "original1",
                                   "MeshKeyOrig2"   : "original2",
                                   "MeshKeyClipped" : "clippedMesh",
                                   "ClippedMarker1"        : 2,
                                   "ClippedMarker2"        : 1,
                                   "Tolerance"      : 1e-8,
                                }),)

dm.DoOps(dataManagerOps) or FSError.PrintAndExit()

#meshOps = ("PrintInfo", "Check")
#fsmeshClipped.DoOps(meshOps) or FSError.PrintAndExit()
fsmeshClipped.PrintInfo()
fsmeshOrig1.PrintInfo()
