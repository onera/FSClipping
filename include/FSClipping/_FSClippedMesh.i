%module FSClipping

%{
// --- specific includes required to compile the wrapper code ---
#include "FSDataManagerOpFactory.h"
#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"

#ifdef FS_WITH_NAMESPACE
using namespace FS;
#endif
%}

// --- import the wrappers we depend on ---
%import(module="FSDataManager") "FSConfig.i"
%import(module="FSDataManager") "FSTypes.i"
%import(module="FSDataManager") "FSString.i"
%import(module="FSDataManager") "FSMeshEnums.i"
%import(module="FSDataManager") "FSDataManagerOpParams.i"

// use contents of SWIG interface file 'std_string.i'
%import(module="FSDataManager") "std_string.i"

// --- gather SWIG interface files ---
%include "FSClippedMeshParams.i"

// --- and now include the stuff for wrapping it in a datamanager-op
%init %{
FSDataManagerOpFactory.Register("ClippedMesh", FSClippedMesh::Create, FSClippedMeshParams::Create);
%}

%pythoncode %{
  """
  Register the mesh operation parameters.
  """
  import sys
  if 'FSDM.FSDataManager' in sys.modules:
    raise ImportError("Do not import FSDataManager as a sub-module of FSDM, plugins are unable to register properly.")
  from FSDataManager import FSDataManagerOpParamsFactory
  FSDataManagerOpParamsFactory.Register("ClippedMesh", FSClippedMeshParams())
%}
