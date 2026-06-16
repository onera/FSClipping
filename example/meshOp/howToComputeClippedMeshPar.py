import FSDM
import os, math
from FSDataManager import FSClac, FSLog, FSError, FSMesh, FSMeshEnums, FSQuantityDescArray, FSUtil, FSFloatArray, FSStringArray
from FSDataManager import FSDataName, FSDataSpecArray, FSDatasetInfo, FSQuantityDesc, FSDataValue
from FSDataManager import FSDataManager

import FSClipping

# --- parallel setup: 3 processes required
# proc 0 = subject mesh (mesh1), proc 1 = clipper mesh (mesh2), proc 2 = matcher
globalClac = FSClac()
procId = globalClac.GetProcID()
meshId = -1
if(procId) :
    if(procId%2):
        meshId = 0
    else:
        meshId = 1
else:
    meshId = 2
if(meshId < 0):
    FSError.PrintAndExit()

print("MeshId " + str(meshId) + " get proc " + str(procId))
clac = FSClac()
globalClac.DivideIntoGroups(meshId, clac)

# --- instantiate data manager with global communicator
dm = FSDataManager(globalClac)

# --- load meshes on their respective processes
if meshId == 0:
    fsmeshOrig1 = dm.GetMesh("original1", clac)
    meshOps = (("ImportMeshTAU", {"MeshFilename": "test/Mesh/input/cube_hexa_fine.grid"}),
               "RepartitionMeshRCB",
               "CreateLocalNumbering",
               "PrintInfo",
               "Check",
              )
    fsmeshOrig1.DoOps(meshOps) or FSError.PrintAndExit()

if meshId == 1:
    fsmeshOrig2 = dm.GetMesh("original2", clac)
    meshOps = (("ImportMeshTAU", {"MeshFilename": "test/Mesh/input/cube_tetra_fine.grid"}),
               "RepartitionMeshRCB",
               "CreateLocalNumbering",
               "PrintInfo",
               "Check",
              )
    fsmeshOrig2.DoOps(meshOps) or FSError.PrintAndExit()


# --- run parallel clipping operation
dataManagerOps = (("ClippedMesh", {"MeshKeyOrig1"      : "original1",
                                   "MeshKeyOrig2"      : "original2",
                                   "MeshKeyClipped1"   : "clippedMesh1",
                                   "MeshKeyClipped2"   : "clippedMesh2",
                                   "ClippedMarker1"    : 2,
                                   "ClippedMarker2"    : 1,
                                   "Tolerance"         : 1e-8,
                                }),)

dm.DoOps(dataManagerOps) or FSError.PrintAndExit()

if meshID == 0:
   if(dm.HasMesh("clippedMesh1")):
       clippedMesh1 = dm.GetMesh("clippedMesh1", False)
       clippedMesh1.PrintInfo()
   else:
      FSError.PrintAndExit()

if meshID == 1:
   if(dm.HasMesh("clippedMesh2")):
       clippedMesh1 = dm.GetMesh("clippedMesh2", False)
       clippedMesh1.PrintInfo()
   else:
      FSError.PrintAndExit()

#if meshId == 0:
#    meshOps = (("ExportMeshVTK", {"Filename"          : "meshHexa",
#                              "Format"            : "RAW",
#                              "FilePerProcess"    : True,
#                              "PrefixDatasetName" : True,
#                              "SplitDataset"      : True,
#                              "ExtractVectors"    : True}),)
#    fsmeshOrig1.DoOps(meshOps)
#
#if meshId == 1:
#    meshOps = (("ExportMeshVTK", {"Filename"          : "meshTetra",
#                              "Format"            : "RAW",
#                              "FilePerProcess"    : True,
#                              "PrefixDatasetName" : True,
#                              "SplitDataset"      : True,
#                              "ExtractVectors"    : True}),)
#    fsmeshOrig2.DoOps(meshOps)



