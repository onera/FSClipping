#ifndef FSCLIPPEDMESH_H
#define FSCLIPPEDMESH_H

// --- include base class ---
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSClippedMeshParams.h"
#include "FSClipping/FSTopologyAssembler.h"
#include "FSDataManagerOp.h"
#include <FSMesh.h>
#include <FSMeshFaceExtractor.h>

// --- include operation parameter class ---

_FS_BEGIN_NAMESPACE

class FSClippedMesh : public FSDataManagerOp
{
  //! Define something typedef

public:
  // --- Class Methods ---

  //! A creation function for a new FSLinearSubElementsMesh-object.
  /*!
    \param clac the parallel communication context for the mesh operation.
    \return A new FSClippedMesh-object.
  */
  static FSDataManagerOp* Create(FSClac* clac);


  // --- Methods ---

  //! A constructor.
  /*!
    \param clac the parallel communication context for the mesh operation.
  */
  FSClippedMesh(FSClac* clac);

  //! A destructor.
  virtual ~FSClippedMesh();

  // ==============================================================
  //
  //                   DataManagerOp Interface
  //

  //! Get the name of the class of the object for the creation of a new object in a factory.
  /*!
    \return The name of the class.
  */
  FSString GetClassName() const override;

  //! Perform the operation on a mesh.
  /*!
    \param data the pointer to the data manager data to be operated on.
    \param params the data manager operation parameters.
    \return True if the operation succeeded, otherwise false.
  */
  bool DoOp(FSDataManagerData*& data, const FSDataManagerOpParams* params) override;

protected:
  // --- Methods ---
  //! A constructor, disabled.
  FSClippedMesh();

  //! A copy constructor, disabled.
  FSClippedMesh(const FSClippedMesh&);
  //! Generate the mesh clipped of the mesh 1 (subject) with the mesh 2 (clipped).
  /*!
    \param meshOriginal1 The original mesh1 to clip.
    \param meshOriginal2 The original mesh2 used for the clipping.
    \param[out] meshClipped The destination mesh of linear sub-elements.
    \return True if the operation succeeded, otherwise false.
  */
  bool GenerateClippedMesh(FSMesh& meshOriginal1, FSMesh& meshOriginal2, FSMesh& clippedMesh);

  //! Extract all the boundary data of mesh 1 and mesh 2 used for the clipping algorithm.
  /*!
    \param meshOriginal1 The original mesh1 to clip.
    \param meshOriginal2 The original mesh2 used for the clipping.
    \param[out] boundaryExtraction1 .
    \param[out] boundaryExtraction2 .
    \return True if the operation succeeded, otherwise false.
  */
  bool ExtractBoundaryFaces(FSMesh& meshOriginal1, FSMesh& meshOriginal2, BoundaryExtraction& boundaryExtraction1, BoundaryExtraction& boundaryExtraction2);

  //! Clipping algorithm.
  /*!
    \param meshOriginal2 The original mesh2 used for the clipping.
    \param boundaryExtraction1 .
    \param boundaryExtraction2 .
    \param[inout] topologyAssembler .
    \param[out] clippedSurface .
    \return True if the operation succeeded, otherwise false.
  */
  bool ComputeSurfaceClipped(FSMesh& meshOriginal2, BoundaryExtraction& boundaryExtraction1, BoundaryExtraction& boundaryExtraction2, FSTopologyAssembler& topologyAssembler, FSTopologyData& clippedSurface);

  //! Generathe the new topology of the mesh 1 with the new polygon faces.
  /*!
    \param meshOriginal1 The original mesh1
    \param boundaryExtraction1 .
    \param topologyAssembler .
    \param surfaceClippedTopo .
    \param[out] meshClippedTopo .
    \return True if the operation succeeded, otherwise false.
  */
  bool GenerateMeshClippedTopo(FSMesh& meshOriginal1, FSMesh& meshOriginal2, BoundaryExtraction& boundaryExtraction1, BoundaryExtraction& boundaryExtraction2, FSTopologyData& meshClippedTopo);

  //! Generate the final clipped mesh.
  /*!
    \param meshOriginal1 The original mesh1
    \param meshClippedTopo .
    \param[out] clippedMesh .
    \return True if the operation succeeded, otherwise false.
  */
  bool GenerateMesh(FSMesh& meshOriginal1, FSTopologyData& meshClippedTopo, FSMesh& clippedMesh);

  // --- Menbers ---

  //! The clipped mesh operation parameters.
  FSClippedMeshParams mParams;

  //! The face extractor to get the faces connecting elements.
  FSMeshFaceExtractor mFaceExtractor1, mFaceExtractor2;
};

_FS_END_NAMESPACE
#endif
