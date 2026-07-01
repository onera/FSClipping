import FSDM
import os, math
from FSDataManager import FSClac, FSLog, FSError, FSMesh, FSMeshEnums, FSQuantityDescArray, FSUtil, FSFloatArray, FSStringArray
from FSDataManager import FSDataName, FSDataSpecArray, FSDatasetInfo, FSQuantityDesc, FSDataValue
from FSDataManager import FSDataManager

import FSClipping

# --- parallel setup: 3 processes required
# proc 0 = subject mesh (mesh1), proc 1 = clipper mesh (mesh2), proc 2 = matcher
globalClac = FSClac()
assert globalClac.GetNProcs() == 3, "3 MPI processes required"
meshID = globalClac.GetProcID()
clac = FSClac()
globalClac.DivideIntoGroups(meshID, clac)

# --- instantiate data manager with global communicator
dm = FSDataManager(globalClac)

# --- load meshes on their respective processes
if meshID == 0:
    fsmeshOrig1 = dm.GetMesh("original1", clac)
    meshOps = (("ImportMeshTAU", {"MeshFilename": "test/Mesh/input/cube_hexa_coarse.grid"}),
               "CreateLocalNumbering",
               "PrintInfo",
               "Check",
              )
    fsmeshOrig1.DoOps(meshOps) or FSError.PrintAndExit()

if meshID == 1:
    fsmeshOrig2 = dm.GetMesh("original2", clac)
    meshOps = (("ImportMeshTAU", {"MeshFilename": "test/Mesh/input/cube_hexa_fine.grid"}),
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
       meshOps = ( "PrintInfo", ("ExportMeshHDF5", {"MeshFilename" : "test/Mesh/output/cube_hexa_coarse_clipped.h5"}),)
       clippedMesh1.DoOps(meshOps) or FSError.PrintAndExit() 
   else:
      FSError.PrintAndExit()

if meshID == 1:
   if(dm.HasMesh("clippedMesh2")):
       clippedMesh2 = dm.GetMesh("clippedMesh2", False)
       meshOps = ( "PrintInfo", ("ExportMeshHDF5", {"MeshFilename" : "test/Mesh/output/cube_hexa_fine_clipped.h5"}),)
       clippedMesh2.DoOps(meshOps) or FSError.PrintAndExit() 
   else:
      FSError.PrintAndExit()


