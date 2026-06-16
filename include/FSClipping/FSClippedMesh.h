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
  //! Top-level parallel driver: dispatches to proc-specific steps based on meshId.
  /*!
    \param meshId  Global proc rank (0 = subject, 1 = clipper, 2 = matcher).
    \param meshClac Single-proc communicator for this process (from DivideIntoGroups).
    \return True if the operation succeeded, otherwise false.
  */
  bool GenerateClippedMesh(FS_intT meshId, FSClac& meshClac);

  //! Extract boundary faces for a single mesh.
  /*!
    \param mesh    The mesh to extract from.
    \param marker  Boundary marker to select.
    \param[out] be Extracted boundary data.
    \return True if any faces were extracted, otherwise false.
  */
  bool ExtractBoundaryFaces(FSMesh& mesh, FS_intT marker, BoundaryExtraction& be);

  //! Build the new topology of mesh1 from pre-computed face matches.
  /*!
    \param mesh1           The subject mesh (proc 0).
    \param be1             Boundary extraction of mesh1.
    \param matches         Face matches received from the matcher proc.
    \param[out] meshClippedTopo Resulting topology.
    \return True if the operation succeeded, otherwise false.
  */
  bool GenerateMeshClippedTopo(FSMesh& mesh1, const BoundaryExtraction& be1,
                                const std::vector<FSFaceMatch>& matches,
                                FSTopologyData& meshClippedTopo);

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
