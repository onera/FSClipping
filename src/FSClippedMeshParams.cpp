#include "FSClipping/FSClippedMeshParams.h"

_FS_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//
//  Create
//

FSDataManagerOpParams* FSClippedMeshParams::Create()
{
  return new FSClippedMeshParams();
}


//-----------------------------------------------------------------------------
//
//  Copy constructor
//

FSClippedMeshParams::FSClippedMeshParams(const FSClippedMeshParams& other)
{
  Copy(other);
}


//-----------------------------------------------------------------------------
//
//  Assignment operator
//

FSClippedMeshParams& FSClippedMeshParams::operator=(const FSClippedMeshParams& other)
{
  Copy(other);
  return *this;
}

//-----------------------------------------------------------------------------
//
//  Copy
//

void FSClippedMeshParams::Copy(const FSClippedMeshParams& other)
{
  // --- base class ---
  FSDataManagerOpParams::Copy(other);

  // --- self ---
  mMeshKeyOriginal1 = other.mMeshKeyOriginal1;
  mMeshKeyOriginal2 = other.mMeshKeyOriginal2;
  mMeshKeyClipped = other.mMeshKeyClipped;
  mMarker1 = other.mMarker1;
  mMarker2 = other.mMarker2;
  mTol = other.mTol;
}

//-----------------------------------------------------------------------------
//
//  Duplicate
//

void FSClippedMeshParams::Duplicate(const FSDataManagerOpParams& other)
{
  Copy(dynamic_cast<const FSClippedMeshParams&>(other));
}


//-----------------------------------------------------------------------------
//
//  Reset
//

void FSClippedMeshParams::Reset()
{
  *this = FSClippedMeshParams{};
}

//-----------------------------------------------------------------------------
//
//  IsInitialized
//

bool FSClippedMeshParams::IsInitialized() const
{
  if(mMeshKeyOriginal1.IsInitialized() && mMeshKeyOriginal2.IsInitialized() && mMarker1 > -1 && mMarker2 > -1) {
    return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
//
//  GetClassName
//

FSString FSClippedMeshParams::GetClassName() const
{
  return "FSClippedMeshParams";
}

#ifdef FS_HAVE_LIBXML2

//-----------------------------------------------------------------------------
//
//  InitFromXML
//
bool FSClippedMeshParams::InitFromXML(xmlDocPtr doc, xmlNodePtr node)
{
  // --- check ---
  if(doc == NULL) {
    FSError.Set("NULL pointer given for XML document, non-blanked mesh parameters can not be initialized from XML node.");
    return false;
  }
  if(node == NULL) {
    FSError.Set("NULL pointer given for XML node, non-blanked mesh parameters can not be initialized from XML node.");
    return false;
  }

  // --- base class ---
  if(!FSDataManagerOpParams::InitFromXML(doc, node)) {
    return false;
  }
  // --- self ---
  bool marker1Parsed = false;
  bool marker2Parsed = false;
  bool tolParsed = false;
  bool meshKeyOrig1Parsed = false;
  bool meshKeyOrig2Parsed = false;

  xmlNodePtr cur = node->children;

  while(cur) {
    if(cur->type == XML_ELEMENT_NODE) {
      // --- check for 'marker' ---
      if(xmlStrcasecmp(cur->name, (const xmlChar*)"ClippedMarker1") == 0 ||
         xmlStrcasecmp(cur->name, (const xmlChar*)"Clipped marker 1") == 0 ||
         xmlStrcasecmp(cur->name, (const xmlChar*)"Marker 1") == 0) {
        // --- check ---
        if(marker1Parsed) {
          FSString msg = "Error during parsing XML node '";
          msg += FSString::UTF8TOISOLAT1(cur->name);
          msg += "'";
          msg += ", ClippedMarker1  already processed.";
          FSError.Set(msg.c_str());
          return false;
        }
        marker1Parsed = true;
        // --- get int ---
        if(!FSXMLUtil::ParseInt(doc, cur, mMarker1)) {
          return false;
        }
      }
      // --- check for 'marker' ---
      else if(xmlStrcasecmp(cur->name, (const xmlChar*)"ClippedMarker2") == 0 ||
              xmlStrcasecmp(cur->name, (const xmlChar*)"Clipped marker 2") == 0 ||
              xmlStrcasecmp(cur->name, (const xmlChar*)"Marker 2") == 0) {
        // --- check ---
        if(marker2Parsed) {
          FSString msg = "Error during parsing XML node '";
          msg += FSString::UTF8TOISOLAT1(cur->name);
          msg += "'";
          msg += ", ClippedMarker2 already processed.";
          FSError.Set(msg.c_str());
          return false;
        }
        marker2Parsed = true;
        // --- get int ---
        if(!FSXMLUtil::ParseInt(doc, cur, mMarker2)) {
          return false;
        }
      }
      // --- check for 'meshkeyOrig1' ---
      else if(xmlStrcasecmp(cur->name, (const xmlChar*)"MeshKey Original1") == 0 ||
              xmlStrcasecmp(cur->name, (const xmlChar*)"MeshKey Orig1") == 0 ||
              xmlStrcasecmp(cur->name, (const xmlChar*)"MeshKeyOrig1") == 0) {
        // --- check ---
        if(meshKeyOrig1Parsed) {
          FSString msg = "Error during parsing XML node '";
          msg += FSString::UTF8TOISOLAT1(cur->name);
          msg += "'";
          msg += ", MeshKeyOrig1 already processed.";
          FSError.Set(msg.c_str());
          return false;
        }
        meshKeyOrig1Parsed = true;
        // --- get string ---
        xmlChar* nodeValue = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
        mMeshKeyOriginal1 = FSString::UTF8TOISOLAT1(nodeValue);
        mMeshKeyOriginal2.Strip();
        // --- cleanup ---
        xmlFree(nodeValue);
      }
      // --- check for 'meshkeyOrig2' ---
      else if(xmlStrcasecmp(cur->name, (const xmlChar*)"MeshKey Original2") == 0 ||
              xmlStrcasecmp(cur->name, (const xmlChar*)"MeshKey Orig2") == 0 ||
              xmlStrcasecmp(cur->name, (const xmlChar*)"MeshKeyOrig2") == 0) {
        // --- check ---
        if(meshKeyOrig2Parsed) {
          FSString msg = "Error during parsing XML node '";
          msg += FSString::UTF8TOISOLAT1(cur->name);
          msg += "'";
          msg += ", MeshKeyOrig2 already processed.";
          FSError.Set(msg.c_str());
          return false;
        }
        meshKeyOrig2Parsed = true;
        // --- get string ---
        xmlChar* nodeValue = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
        mMeshKeyOriginal2 = FSString::UTF8TOISOLAT1(nodeValue);
        mMeshKeyOriginal2.Strip();
        // --- cleanup ---
        xmlFree(nodeValue);
      }

      // --- check for 'marker' ---
      else if(xmlStrcasecmp(cur->name, (const xmlChar*)"ClippedTolerance") == 0 ||
              xmlStrcasecmp(cur->name, (const xmlChar*)"Clipped tolerance") == 0 ||
              xmlStrcasecmp(cur->name, (const xmlChar*)"Tolerance") == 0) {
        // --- check ---
        if(tolParsed) {
          FSString msg = "Error during parsing XML node '";
          msg += FSString::UTF8TOISOLAT1(cur->name);
          msg += "'";
          msg += ", ClippedMarker2 already processed.";
          FSError.Set(msg.c_str());
          return false;
        }
        tolParsed = true;
        // --- get int ---
        if(!FSXMLUtil::ParseFloat(doc, cur, mTol)) {
          return false;
        }
      } else {
        // ... ignore ...
      }
    }
    cur = cur->next;
  }

  return true;
}

//-----------------------------------------------------------------------------
//
//  ConvertToXML
//
bool FSClippedMeshParams::ConvertToXML(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr node)
{
  // --- check ---
  if(doc == nullptr) {
    FSError.Set("NULL pointer given for XML document, data transformation parameters can not be converted to an XML node.");
    return false;
  }
  if(node == nullptr) {
    FSError.Set("NULL pointer given for XML node, data transformation parameters can not be be converted to an XML node.");
    return false;
  }

  // --- base class ---
  if(!FSDataManagerOpParams::ConvertToXML(doc, ns, node)) {
    return false;
  }

  // --- self ---

  // --- convert mesh keys ---
  if(mMeshKeyOriginal1.IsInitialized()) {
    xmlNewChild(node, ns, reinterpret_cast<const xmlChar*>("MeshKeyOriginal1"), reinterpret_cast<const xmlChar*>(mMeshKeyOriginal1.c_str()));
  }

  if(mMeshKeyOriginal2.IsInitialized()) {
    xmlNewChild(node, ns, reinterpret_cast<const xmlChar*>("MeshKeyOriginal2"), reinterpret_cast<const xmlChar*>(mMeshKeyOriginal2.c_str()));
  }
  return true;
}
#endif

//-----------------------------------------------------------------------------
//
//  GetBufSize
//

FSClippedMeshParams::sizeT FSClippedMeshParams::GetBufSize(FSClac& clac) const
{
  sizeT bufSize = 0;

  // --- base class ---
  bufSize += FSDataManagerOpParams::GetBufSize(clac);

  // --- self ---
  bufSize += clac.GetBufSize<FS_intT>();   // mMarker1
  bufSize += clac.GetBufSize<FS_intT>();   // mMarker2
  bufSize += clac.GetBufSize<FS_floatT>(); // /Tol
  bufSize += mMeshKeyOriginal1.GetBufSize(clac);
  bufSize += mMeshKeyOriginal2.GetBufSize(clac);
  bufSize += mMeshKeyClipped.GetBufSize(clac);

  return bufSize;
}

//-----------------------------------------------------------------------------
//
//  Pack
//

void FSClippedMeshParams::Pack(FSClac& clac)
{
  // --- base class ---
  FSDataManagerOpParams::Pack(clac);

  // --- self ---
  clac.Pack(&mMarker1);
  clac.Pack(&mMarker2);
  clac.Pack(&mTol);
  mMeshKeyOriginal1.Pack(clac);
  mMeshKeyOriginal2.Pack(clac);
  mMeshKeyClipped.Pack(clac);
}


//-----------------------------------------------------------------------------
//
//  Unpack
//

void FSClippedMeshParams::Unpack(FSClac& clac)
{
  // --- base class ---
  FSDataManagerOpParams::Unpack(clac);

  // --- self ---
  clac.Unpack(&mMarker1);
  clac.Unpack(&mMarker2);
  clac.Unpack(&mTol);
  mMeshKeyOriginal1.Unpack(clac);
  mMeshKeyOriginal2.Unpack(clac);
  mMeshKeyClipped.Unpack(clac);
}



_FS_END_NAMESPACE