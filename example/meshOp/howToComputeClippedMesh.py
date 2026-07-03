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
    meshOps = (("ImportMeshTAU", {"MeshFilename": "test/Mesh/input/mesh_cylinder_1.grid"}),
               "CreateLocalNumbering",
               "PrintInfo",
               "Check",
              )
    fsmeshOrig1.DoOps(meshOps) or FSError.PrintAndExit()

if meshID == 1:
    fsmeshOrig2 = dm.GetMesh("original2", clac)
    meshOps = (("ImportMeshTAU", {"MeshFilename": "test/Mesh/input/mesh_cylinder_2.grid"}),
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

if meshID == 0:
   if(dm.HasMesh("clippedMesh1")):
       clippedMesh1 = dm.GetMesh("clippedMesh1", False)
       meshOps = ( "PrintInfo", ("ExportMeshHDF5", {"MeshFilename" : "test/Mesh/output/mesh_cylinder_1_clipped.h5"}),)
       clippedMesh1.DoOps(meshOps) or FSError.PrintAndExit() 
   else:
      FSError.PrintAndExit()

if meshID == 1:
   if(dm.HasMesh("clippedMesh2")):
       clippedMesh2 = dm.GetMesh("clippedMesh2", False)
       meshOps = ( "PrintInfo", ("ExportMeshHDF5", {"MeshFilename" : "test/Mesh/output/mesh_cylinder_2_clipped.h5"}),)
       clippedMesh2.DoOps(meshOps) or FSError.PrintAndExit() 
   else:
      FSError.PrintAndExit()

if meshID == 0:
    meshOps = (("ExportMeshVTK", {"Filename"            : "test/Mesh/output/mesh_clipped_cylinder_1",
                              "Format"            : "RAW",
                              "FilePerProcess"    : True,
                              "PrefixDatasetName" : True,
                              "SplitDataset"      : True,
                              "ExtractVectors"    : True}),)

    if not clippedMesh1.DoOps(meshOps):
        FSError.PrintAndExit()
        
if meshID == 1:
    meshOps = (("ExportMeshVTK", {"Filename"          : "test/Mesh/output/mesh_clipped_cylinder_2",
                              "Format"            : "RAW",
                              "FilePerProcess"    : True,
                              "PrefixDatasetName" : True,
                              "SplitDataset"      : True,
                              "ExtractVectors"    : True}),)

    if not clippedMesh2.DoOps(meshOps):
        FSError.PrintAndExit()