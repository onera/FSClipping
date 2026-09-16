import FSDM
import os, math
from FSDataManager import FSClac, FSLog, FSError, FSMesh, FSMeshEnums, FSQuantityDescArray, FSUtil, FSFloatArray, FSStringArray
from FSDataManager import FSDataName, FSDataSpecArray, FSDatasetInfo, FSQuantityDesc, FSDataValue
from FSDataManager import FSDataManager

import FSClipping

def MeshID(procID, nMesh1):
    if(procID == 0):
        return 0
    elif(procID <= nMesh1):
        return 1
    return 2

# --- parallel setup
globalClac = FSClac()
nTotalProc = globalClac.NWorldProcs()
nProcMesh1 = 3
nProcMesh2 = 2

if(nProcMesh1 + nProcMesh2 != nTotalProc - 1) :
    FSError.PrintAndExit("Wrong number of mpi process were chosen.")

procID = globalClac.WorldProcID()
meshID = MeshID(procID, nProcMesh1)

clac = FSClac()
globalClac.DivideIntoGroups(meshID, clac)

# --- instantiate data manager with global communicator
dm = FSDataManager(globalClac)

# --- load meshes on their respective processes
if meshID == 1:
    fsmeshOrig1 = dm.GetMesh("original1", clac)
    meshOps = (("ImportMeshTAU", {"MeshFilename": "test/Mesh/input/mesh_cylinder_1.grid"}),
               "RepartitionMeshRCB",
               "CreateLocalNumbering",
               "PrintInfo",
               "Check",
              )
    fsmeshOrig1.DoOps(meshOps) or FSError.PrintAndExit()

if meshID == 2:
    fsmeshOrig2 = dm.GetMesh("original2", clac)
    meshOps = (("ImportMeshTAU", {"MeshFilename": "test/Mesh/input/mesh_cylinder_2.grid"}),
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
                                   "ClippedMarker1"    : 1,
                                   "ClippedMarker2"    : 1,
                                   "Tolerance"         : 1e-8,
                                }),)

dm.DoOps(dataManagerOps) or FSError.PrintAndExit()

if meshID == 1:
   if(dm.HasMesh("clippedMesh1")):
       clippedMesh1 = dm.GetMesh("clippedMesh1", False)
       meshOps = ( "PrintInfo",)
       clippedMesh1.DoOps(meshOps) or FSError.PrintAndExit() 
   else:
      FSError.PrintAndExit()

if meshID == 2:
   if(dm.HasMesh("clippedMesh2")):
       clippedMesh2 = dm.GetMesh("clippedMesh2", False)
       meshOps = ( "PrintInfo",)
       clippedMesh2.DoOps(meshOps) or FSError.PrintAndExit() 
   else:
      FSError.PrintAndExit()