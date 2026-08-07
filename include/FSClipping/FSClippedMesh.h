#ifndef FSCLIPPEDMESH_H
#define FSCLIPPEDMESH_H

#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSClippedMeshParams.h"
#include "FSDataManagerOp.h"
#include <FSMesh.h>
#include <FSMeshFaceExtractor.h>

_FS_BEGIN_NAMESPACE

class FSClippedMesh : public FSDataManagerOp
{
public:
  //! Factory function — returns a new FSClippedMesh instance.
  static FSDataManagerOp* Create(FSClac* clac);

  //! Constructor.
  FSClippedMesh(FSClac* clac);

  //! Destructor.
  virtual ~FSClippedMesh();

  // ==============================================================
  //                   DataManagerOp interface
  // ==============================================================

  FSString GetClassName() const override;

  //! Runs the clipping operation on the two meshes referenced by params.
  bool DoOp(FSDataManagerData*& data, const FSDataManagerOpParams* params) override;

protected:
  FSClippedMesh();
  FSClippedMesh(const FSClippedMesh&);

  //! Top-level parallel driver: dispatches to proc-specific steps.
  //! Proc 0 = subject, proc 1 = clipper, proc 2 = matcher.
  bool GenerateClippedMesh(FS_intT meshId);

  // --- Members ---

  FSClippedMeshParams mParams;
};

_FS_END_NAMESPACE
#endif
