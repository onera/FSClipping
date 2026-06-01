%{
// --- includes required to compile the wrapper code ---
#include "FSClipping/FSClippedMeshParams.h"
%}

// --- define declarations to be ignored ---
%ignore FSClippedMeshParams::operator=(const FSClippedMeshParams&);
%ignore FSClippedMeshParams::GetBufSize;
%ignore FSClippedMeshParams::Pack;
%ignore FSClippedMeshParams::Unpack;

// --- define ANSI C/C++ declarations to be interfaced ---
%include "FSClippedMeshParams.h"

/* extend FSClippedMeshParams for a more easy usage */
%extend FSClippedMeshParams {
  %pythoncode %{
  def ToXML(self, params):
      """
      Convert data transformation parameters defined by a nested data structure
      into an XML representation.
      """
      if type(params) != dict:
          raise Exception("not a mapping type, %s can not be converted into XML format for class %s!" % (str(params), str(__class__)))
      # --- base-class ---
      paramsXML = super(FSClippedMeshParams, self).ToXML(params)

      # --- self ---
      for key in ("MeshKey Original 1",
                  "MeshKey Orig 1",
                  "MeshKeyOrig1"):
          if key in params:
              paramsXML += "<%s>%s</%s>\n" % (key, str(params[key]), key)
              del params[key]

      for key in ("MeshKey Original 2",
                  "MeshKey Orig 2",
                  "MeshKeyOrig2"):
          if key in params:
              paramsXML += "<%s>%s</%s>\n" % (key, str(params[key]), key)
              del params[key]

      for key in ("MeshKeyClipped1",):
          if key in params:
              paramsXML += "<%s>%s</%s>\n" % (key, str(params[key]), key)
              del params[key]

      for key in ("MeshKeyClipped2",):
          if key in params:
              paramsXML += "<%s>%s</%s>\n" % (key, str(params[key]), key)
              del params[key]

      for key in ("Clipped marker 1",
                  "ClippedMarker1",
                  "Marker1"):
          if key in params:
              paramsXML += "<%s>%s</%s>\n" % (key, str(params[key]), key)
              del params[key]

      for key in ("Clipped marker 2",
                  "ClippedMarker2",
                  "Marker2"):
          if key in params:
              paramsXML += "<%s>%s</%s>\n" % (key, str(params[key]), key)
              del params[key]

      for key in ("Clipped tolerance",
                  "ClippedTolerance",
                  "Tolerance"):
          if key in params:
              paramsXML += "<%s>%s</%s>\n" % (key, str(params[key]), key)
              del params[key]

      return paramsXML
  %}
};
