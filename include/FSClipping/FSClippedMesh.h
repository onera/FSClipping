#ifndef FSCLIPPEDMESH_H
#define FSCLIPPEDMESH_H

#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSClippedMeshParams.h"
#include "FSClipping/FSTopologyAssembler.h"
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
  bool GenerateClippedMesh(FS_intT meshId, FSClac& meshClac);

  //! Extract boundary faces for the given mesh and boundary marker.
  bool ExtractBoundaryFaces(FSMesh& mesh, FS_intT marker, BoundaryExtraction& be);

  //! Build the clipped topology of mesh1 from pre-computed face matches.
  bool GenerateMeshClippedTopo(FSMesh& mesh1, const BoundaryExtraction& be1, const std::vector<FSFaceMatch>& matches,
                               FSTopologyData& meshClippedTopo);

  //! Reconstruct the final clipped mesh from the topology data.
  bool GenerateMesh(FSMesh& meshOriginal1, FSTopologyData& meshClippedTopo, FSMesh& clippedMesh);

  // --- Members ---

  FSClippedMeshParams mParams;
  FSMeshFaceExtractor mFaceExtractor1, mFaceExtractor2;
};

_FS_END_NAMESPACE
#endif
