#include "FSClipping/FSClippingEngine.h"
#include "FSClipping/FSCellClipped.h"
#include "FSClipping/FSCellSource.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSMeshReconstruction.h"
#include "FSClipping/FSNodeRegistry.h"
#include "FSClipping/FSTopologySink.h"
#include <FSMeshData.h>
#include <numeric>

_FS_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//
//  BuildTopology
//

FSTopologyData FSClippingEngine::BuildTopology(FSMesh& mesh, const BoundaryExtraction& be,
                                               const std::vector<FSFaceMatch>& matches,
                                               const std::vector<FSNodeMove>& moves)
{
  // 1 - Build the input cells
  const std::vector<FSCellSource> sources = BuildCellSources(mesh);
  const FSCellClipped clipped(matches, be);

  // 2 - Build the nodes
  const FSNodeRegistry nodes(mesh, sources, clipped, tol_, moves);

  // 3 - Build the topology
  FSTopologySink out(nodes);

  for(const auto& src : sources) // Cell source first
    src.ForEachOwnedCell([&](FS_intT c) {
      if(clipped.HasBoundaryFace(c))
        out.AddPoly3D(c, src.Type(), BuildClippedCell(src, c, clipped, nodes));
      else if(clipped.IsBoundarySurface(c))
        return;
      else if(src.IsPoly3D())
        out.AddPoly3D(src, c);
      else if(src.IsPoly2D()) // marker -1: read back from the cell itself
        out.AddPoly2D(c, src.Type(), src.Row(c), -1);
      else
        out.AddStandard(src, c);
    });

  clipped.ForEachMatch([&](FS_intT m, const FSFaceMatch& f) {
    out.AddPoly2D(f.faceOwner1, f.faceOwnerType1, nodes.Polygon(m), f.faceMarker1);
  });

  return std::move(out).Finalize();
}

//-----------------------------------------------------------------------------
//
//  Reconstruct
//

FSMesh FSClippingEngine::Reconstruct(FSClac& clac, FSMesh& meshOriginal, FSTopologyData& topo, bool copyAttributes)
{
  FSMeshReconstruction meshReconstruction(clac);
  FSMesh clippedMesh = meshReconstruction.Build(meshOriginal.GetMeshData()->GetUnstructCells(), topo);
  if(copyAttributes)
    meshReconstruction.CopyAttributes(meshOriginal, clippedMesh);
  return clippedMesh;
}

//-----------------------------------------------------------------------------
//
//  RunMatcherProc
//

bool FSClippingEngine::RunMatcherProc(FSClac& globalClac, FSClac& localClac, FS_floatT tol, FS_floatT snapTol)
{
  auto allGathered = FSFaceExchange::GatherReceiveAll(globalClac);
  if(allGathered.meshA.faces.faces.empty() || allGathered.meshB.faces.faces.empty())
    return false;

  // Mesh B interface nodes nearly coincident with mesh A ones are moved onto them before the clipping
  const std::vector<FSNodeMove> moves =
    FSNodeSnap::Compute(allGathered.meshA.faces.faces, allGathered.meshB.faces.faces, snapTol);
  FSNodeSnap::Apply(moves, allGathered.meshB.faces.faces);
  if(!moves.empty()) {
    FSString msg = "FSClippingEngine: ";
    msg.Add(FS_intT(moves.size()));
    msg += " interface nodes of mesh 2 snapped onto mesh 1\n";
    FSLog(&localClac, 0, msg.c_str(), 0);
  }

  FSFaceMatcher matcher(localClac, allGathered.meshA.faces.faces, allGathered.meshB.faces.faces, tol);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  if(matches.empty()) {
    return false;
  }

  // Each side gets its own global numbering (nodes + Poly2D cells)
  FSFaceMatcher::AssignGlobalNodeIds(matches, tol);
  FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshA);

  matcher.ComputeInvertedMatches(matches);
  FSFaceMatcher::AssignGlobalNodeIds(matches, tol);
  FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshB, moves);
  return true;
}

_FS_END_NAMESPACE
