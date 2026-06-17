#ifndef FSCLIPPEDMESHPARAMES_H
#define FSCLIPPEDMESHPARAMES_H

#include "FSDataManagerOpParams.h"

_FS_BEGIN_NAMESPACE

//! Parameters for the FSClippedMesh operation.
class FSClippedMeshParams : public FSDataManagerOpParams
{
public:
  // --- Public parameters ---

  //! Boundary marker of mesh 1 (selects faces to clip).
  FS_intT mMarker1;

  //! Boundary marker of mesh 2 (selects faces to clip).
  FS_intT mMarker2;

  //! Geometric tolerance for the clipping algorithm.
  FS_floatT mTol;

  //! Mesh key of the original mesh 1.
  FSString mMeshKeyOriginal1;

  //! Mesh key of the original mesh 2.
  FSString mMeshKeyOriginal2;

  //! Mesh key for the resulting clipped mesh 1.
  FSString mMeshKeyClipped1;

  //! Mesh key for the resulting clipped mesh 2.
  FSString mMeshKeyClipped2;

  // --- Factory ---

  //! Returns a new FSClippedMeshParams instance.
  static FSDataManagerOpParams* Create();

  // --- Lifecycle ---

  FSClippedMeshParams() = default;
  FSClippedMeshParams(const FSClippedMeshParams&);
  FSClippedMeshParams& operator=(const FSClippedMeshParams& other);

  void Reset() override;
  void Duplicate(const FSDataManagerOpParams& other) override;

  // ==============================================================
  //                        Queries
  // ==============================================================

  bool IsInitialized() const override;
  FSString GetClassName() const override;

  // ==============================================================
  //                          XML
  // ==============================================================

#ifdef FS_HAVE_LIBXML2
  bool InitFromXML(xmlDocPtr doc, xmlNodePtr node) override;
  bool ConvertToXML(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr node) override;
#endif

  // ==============================================================
  //                   Parallel transmission
  // ==============================================================

  sizeT GetBufSize(FSClac& clac) const override;
  void Pack(FSClac& clac) override;
  void Unpack(FSClac& clac) override;

protected:
  void Copy(const FSClippedMeshParams& other);
};

_FS_END_NAMESPACE

#endif // FSCLIPPEDMESHPARAMS_H
