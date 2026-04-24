# --- from the FSDM data manager import the relevant classes you need ---
from FSDataManager import FSDataManagerOpParamsFactory 
from FSDataManager import FSDataManagerOpFactory
# --- import our PlugIn stuff ---
import FSClipping

print("FSClippedMeshParams is registered : ", FSDataManagerOpParamsFactory.IsRegistered("ClippedMesh")) 
print("FSClippedMesh is registered :", FSDataManagerOpFactory.IsRegistered("ClippedMesh"))

