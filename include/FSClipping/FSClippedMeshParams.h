#ifndef FSCLIPPEDMESHPARAMES_H
#define FSCLIPPEDMESHPARAMES_H

#include "FSDataManagerOpParams.h"


_FS_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//
//  FSClippedMeshParams
//

//! This class administrates parameters for the creation of a linear sub-elements mesh.
class FSClippedMeshParams
  : public FSDataManagerOpParams
{
public:
  // --- Member variables (public) ---

  //! \name OCL Parameters
  //@{

  //! The boundary marker of the mesh 1 for the resulting boundary faces in extraction mode.
  FS_intT mMarker1;

  //! The boundary marker of the mesh 2 for the resulting boundary faces in extraction mode.
  FS_intT mMarker2;

  //! The tolerance of the clipping algorithm.
  FS_floatT mTol;

  //! The mesh key of the original mesh 1.
  FSString mMeshKeyOriginal1;

  //! The mesh key of the original mesh 2.
  FSString mMeshKeyOriginal2;

  //! The mesh key of the clipped mesh1.
  FSString mMeshKeyClipped1;

  //! The mesh key of the clipped mesh2.
  FSString mMeshKeyClipped2;
  //@}

  // --- Class Methods ---

  //! A creation function for a new FSClippedMeshParams-object.
  /*!
    \return A new FSClippedMeshParams-object.
  */
  static FSDataManagerOpParams* Create();

  // --- Methods ---

  //! A constructor.
  FSClippedMeshParams() = default;

  //! A copy constructor.
  FSClippedMeshParams(const FSClippedMeshParams&);

  //! An assignment operator.
  FSClippedMeshParams& operator=(const FSClippedMeshParams& other);

  //! Reset all data.
  void Reset() override;

  //! Duplicate the contents of the object.
  /*!
    \param other the instance to be duplicated.
  */
  void Duplicate(const FSDataManagerOpParams& other) override;


  // ==============================================================
  //
  //                       Queries
  //

  //! Check if the parameters are initialized.
  /*!
    \return True if there are any items defined, otherwise false.
  */
  bool IsInitialized() const override;

  //! Get the name of the class of the object for the creation of a new object in a factory.
  /*!
    \return The name of the class.
  */
  FSString GetClassName() const override;


  // ==============================================================
  //
  //                          XML
  //

#ifdef FS_HAVE_LIBXML2
  //   //! Init the entity from an XML node and its children.
  //   /*!
  //     \param doc the current XML document in which the node is included.
  //     \param node the current XML node to be used to init the object.
  //     \return True if the object could be initialized successfully, otherwise false.
  //   */
  bool InitFromXML(xmlDocPtr doc, xmlNodePtr node) override;
  //
  //   //! Convert the entity to an XML node.
  //   /*!
  //     \param doc the current XML document.
  //     \param ns the current XML namespace to be used.
  //     \param node the current XML node to be used to convert the object.
  //     \return True if the object could be converted successfully, otherwise false.
  //   */
  bool ConvertToXML(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr node) override;
#endif


  // ==============================================================
  //
  //                   Parallel Transmission
  //

  //! Get the number of bytes required to pack the object into a buffer.
  /*!
    \param clac the parallel communication context to be used to calculate the buffer size.
    \return The number of bytes required to pack the object into a buffer.
  */
  sizeT GetBufSize(FSClac& clac) const override;

  //! Pack the object into a buffer.
  /*!
    \param clac the parallel communication context which contains the buffer into which the object is packed.
  */
  void Pack(FSClac& clac) override;

  //! Unpack the object from a buffer.
  /*!
    \param clac the parallel communication context which contains the buffer from which the object is unpacked.
  */
  void Unpack(FSClac& clac) override;

protected:
  // --- Methods ---

  //! Perform a deep copy of the object.
  /*!
    \param other the instance to be copied.
  */
  void Copy(const FSClippedMeshParams& other);
};

_FS_END_NAMESPACE

#endif // FSCLIPPEDMESHPARAMS_H