#include "FSClipping/FSVertexSnap.h"
#include "FSClipping/FSFaceMatcher.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <unordered_set>

_FS_BEGIN_NAMESPACE

namespace {

constexpr FS_floatT kInf = std::numeric_limits<FS_floatT>::infinity();

using CellKey = std::array<std::int64_t, 3>;

struct CellKeyHash {
  std::size_t operator()(const CellKey& k) const noexcept
  {
    std::size_t seed = 0;
    for(const auto v : k)
      seed ^= std::hash<std::int64_t>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

// Grid of step snapTol: two points closer than snapTol are in neighbouring cells
CellKey CellOf(const FSVec3& p, FS_floatT step)
{
  return {static_cast<std::int64_t>(std::floor(p[0] / step)), static_cast<std::int64_t>(std::floor(p[1] / step)),
          static_cast<std::int64_t>(std::floor(p[2] / step))};
}

FS_floatT DistanceToSegment(const FSVec3& p, const FSVec3& a, const FSVec3& b)
{
  const FSVec3 ab = b - a;
  const FS_floatT len2 = ab.InnerProduct(ab);
  const FS_floatT t = len2 > 0.0 ? std::clamp((p - a).InnerProduct(ab) / len2, 0.0, 1.0) : 0.0;
  const FSVec3 q(a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]);
  return q.Distance(p);
}

void SortUnique(std::vector<FS_intT>& v)
{
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
}

void Append(std::vector<FS_intT>& to, const std::vector<FS_intT>& from) { to.insert(to.end(), from.begin(), from.end()); }

FS_intT Root(std::vector<FS_intT>& parent, FS_intT i)
{
  while(parent[i] != i) {
    parent[i] = parent[parent[i]];
    i = parent[i];
  }
  return i;
}

std::array<std::uint64_t, 3> KeyOf(const FSVec3& p)
{
  auto bits = [](FS_floatT x) { return std::bit_cast<std::uint64_t>(x == 0.0 ? 0.0 : x); };
  return {bits(p[0]), bits(p[1]), bits(p[2])};
}

void Fail(const char* what, const FSVec3& p)
{
  FSString msg = "FSVertexSnap: ";
  msg += what;
  msg += " near (";
  msg.Add(p[0]);
  msg += ", ";
  msg.Add(p[1]);
  msg += ", ";
  msg.Add(p[2]);
  msg += ")";
  FSError.SetAndPrintAndExit(msg.c_str());
}

} // namespace

std::size_t FSVertexSnap::ExactKeyHash::operator()(const ExactKey& k) const noexcept
{
  std::size_t seed = 0;
  for(const auto v : k)
    seed ^= std::hash<std::uint64_t>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

FS_intT FSVertexSnap::Find(const FSVec3& p) const
{
  const auto it = indexOf_.find(KeyOf(p));
  if(it == indexOf_.end())
    Fail("vertex missing from the raw matches", p);
  return it->second;
}

FSVertexSnap::FSVertexSnap(const std::vector<FSClippingFace>& faces1, const std::vector<FSClippingFace>& faces2,
                           const std::vector<FSFaceMatch>& rawMatches, FS_floatT tol, FS_floatT snapTol)
{
  // 1 - Points, with per side the faces they belong to and the face edges they lie on (within tol)
  struct PointInfo {
    std::array<bool, 2> original{false, false};
    std::array<FS_floatT, 2> outlineDist{kInf, kInf};
    std::array<std::vector<FS_intT>, 2> faces, edges;
  };
  std::vector<PointInfo> info;
  auto addPoint = [&](const FSVec3& p) {
    const auto [it, inserted] = indexOf_.try_emplace(KeyOf(p), FS_intT(points_.size()));
    if(inserted) {
      points_.push_back(p);
      info.emplace_back();
    }
    return it->second;
  };

  std::map<std::pair<ExactKey, ExactKey>, FS_intT> edgeIds;
  auto addMembership = [&](FS_intT s, FS_intT i, const FSClippingFace& f) {
    PointInfo& pi = info[i];
    pi.faces[s].push_back(f.faceIndex());
    const auto& v = f.vertices();
    for(std::size_t k = 0; k < v.size(); ++k) {
      const FS_floatT d = DistanceToSegment(points_[i], v[k], v[(k + 1) % v.size()]);
      pi.outlineDist[s] = std::min(pi.outlineDist[s], d);
      if(d > tol)
        continue;
      ExactKey ka = KeyOf(v[k]), kb = KeyOf(v[(k + 1) % v.size()]);
      if(kb < ka)
        std::swap(ka, kb);
      pi.edges[s].push_back(edgeIds.try_emplace({ka, kb}, FS_intT(edgeIds.size())).first->second);
    }
  };

  const std::array<const std::vector<FSClippingFace>*, 2> faces{&faces1, &faces2};
  std::array<std::unordered_map<FS_intT, const FSClippingFace*>, 2> faceOf;
  for(FS_intT s = 0; s < 2; ++s)
    for(const auto& f : *faces[s]) {
      faceOf[s].emplace(f.faceIndex(), &f);
      for(const auto& v : f.vertices()) {
        const FS_intT i = addPoint(v);
        info[i].original[s] = true;
        addMembership(s, i, f);
      }
    }
  auto addPiece = [&](FS_intT s, FS_intT face, const std::vector<FSVec3>& poly) {
    const auto it = faceOf[s].find(face);
    if(it == faceOf[s].end())
      return;
    for(const auto& p : poly)
      addMembership(s, addPoint(p), *it->second);
  };
  for(const auto& m : rawMatches) {
    addPiece(0, m.face1, m.clippedPoly3D);
    addPiece(1, m.face2, m.clippedPoly3D_face2);
  }
  const FS_intT n = FS_intT(points_.size());

  // 2 - Pairs closer than snapTol, by increasing distance
  struct Pair {
    FS_floatT d;
    FS_intT i, j;
  };
  std::vector<Pair> pairs;
  std::unordered_map<CellKey, std::vector<FS_intT>, CellKeyHash> grid;
  for(FS_intT i = 0; i < n; ++i) {
    const CellKey c = CellOf(points_[i], snapTol);
    for(std::int64_t dx = -1; dx <= 1; ++dx)
      for(std::int64_t dy = -1; dy <= 1; ++dy)
        for(std::int64_t dz = -1; dz <= 1; ++dz) {
          const auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
          if(it == grid.end())
            continue;
          for(const FS_intT j : it->second) {
            const FS_floatT d = points_[i].Distance(points_[j]);
            if(d <= snapTol)
              pairs.push_back({d, i, j});
          }
        }
    grid[c].push_back(i);
  }
  std::stable_sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) { return a.d < b.d; });

  // 3 - Points closer than tol are the same vertex (as for AssignGlobalNodeIds): what they touch is pooled
  std::vector<FS_intT> parent(n);
  std::iota(parent.begin(), parent.end(), 0);
  for(const Pair& p : pairs) {
    if(p.d > tol)
      break;
    parent[Root(parent, p.i)] = Root(parent, p.j);
  }
  std::vector<FS_intT> vertexOf(n);
  std::unordered_map<FS_intT, FS_intT> vertexOfRoot;
  for(FS_intT i = 0; i < n; ++i) {
    const auto [it, inserted] = vertexOfRoot.try_emplace(Root(parent, i), FS_intT(vertices_.size()));
    if(inserted)
      vertices_.emplace_back().outlineDist = {kInf, kInf};
    FSSnapVertex& v = vertices_[it->second];
    vertexOf[i] = it->second;
    const bool anyOriginal = info[i].original[0] || info[i].original[1];
    if(v.point < 0 || anyOriginal)
      v.point = i;
    for(FS_intT s = 0; s < 2; ++s) {
      v.original[s] = v.original[s] || info[i].original[s];
      v.outlineDist[s] = std::min(v.outlineDist[s], info[i].outlineDist[s]);
      Append(v.faces[s], info[i].faces[s]);
      Append(v.edges[s], info[i].edges[s]);
    }
  }
  for(auto& v : vertices_)
    for(FS_intT s = 0; s < 2; ++s) {
      SortUnique(v.faces[s]);
      SortUnique(v.edges[s]);
    }

  // 4 - Greedy merge of the vertices, nearest first; refused when the group would exceed snapTol in diameter or
  // have no valid representative on one side
  const FS_intT nv = FS_intT(vertices_.size());
  std::vector<FS_intT> vparent(nv);
  std::iota(vparent.begin(), vparent.end(), 0);
  std::vector<std::vector<FS_intT>> members(nv);
  for(FS_intT v = 0; v < nv; ++v)
    members[v] = {v};

  std::vector<FS_intT> merged;
  for(const Pair& p : pairs) {
    if(p.d <= tol)
      continue;
    FS_intT a = Root(vparent, vertexOf[p.i]), b = Root(vparent, vertexOf[p.j]);
    if(a == b)
      continue;
    bool tooWide = false;
    for(const FS_intT u : members[a]) {
      for(const FS_intT w : members[b])
        if(points_[vertices_[u].point].Distance(points_[vertices_[w].point]) > snapTol) {
          tooWide = true;
          break;
        }
      if(tooWide)
        break;
    }
    merged = members[a];
    Append(merged, members[b]);
    if(tooWide || Representative(0, merged) == kNoRep || Representative(1, merged) == kNoRep) {
      ++nRefused_;
      continue;
    }
    if(members[a].size() < members[b].size())
      std::swap(a, b);
    vparent[b] = a;
    members[a] = std::move(merged);
    members[b].clear();
  }

  // 5 - Groups and their representative point on each side
  std::unordered_map<FS_intT, FS_intT> groupOfRoot;
  for(FS_intT v = 0; v < nv; ++v)
    groupOfRoot.try_emplace(Root(vparent, v), FS_intT(groupOfRoot.size()));
  group_.resize(n);
  for(FS_intT i = 0; i < n; ++i)
    group_[i] = groupOfRoot.at(Root(vparent, vertexOf[i]));

  for(FS_intT s = 0; s < 2; ++s) {
    rep_[s].assign(groupOfRoot.size(), -1);
    for(const auto& [root, g] : groupOfRoot) {
      const FS_intT r = Representative(s, members[root]);
      rep_[s][g] = r >= 0 ? vertices_[r].point : -1;
    }
  }
}

FS_intT FSVertexSnap::Representative(FS_intT s, const std::vector<FS_intT>& group) const
{
  // The representative lies on every face and face edge the group touches: a point then only slides along its
  // own edge or inside its own face, so the pieces keep the outline of their original face and leave no hole
  std::vector<FS_intT> facesAll, edgesAll;
  FS_intT original = -1;
  bool onSide = false;
  for(const FS_intT v : group) {
    Append(facesAll, vertices_[v].faces[s]);
    Append(edgesAll, vertices_[v].edges[s]);
    onSide = onSide || !vertices_[v].faces[s].empty();
    if(vertices_[v].original[s]) {
      if(original >= 0) // two nodes of one mesh: an edge collapse
        return kNoRep;
      original = v;
    }
  }
  if(!onSide)
    return kNotOnSide;
  SortUnique(facesAll);
  SortUnique(edgesAll);

  auto valid = [&](FS_intT v) {
    const auto& f = vertices_[v].faces[s];
    const auto& e = vertices_[v].edges[s];
    return std::includes(f.begin(), f.end(), facesAll.begin(), facesAll.end()) &&
           std::includes(e.begin(), e.end(), edgesAll.begin(), edgesAll.end());
  };
  if(original >= 0) // an original node never moves
    return valid(original) ? original : kNoRep;

  FS_intT best = kNoRep;
  for(const FS_intT v : group)
    if(!vertices_[v].faces[s].empty() && valid(v) &&
       (best < 0 || vertices_[v].outlineDist[s] < vertices_[best].outlineDist[s]))
      best = v;
  return best;
}

FS_intT FSVertexSnap::Apply(std::vector<FSFaceMatch>& matches, FS_intT side) const
{
  FS_intT dropped = 0;
  std::vector<FSFaceMatch> kept;
  kept.reserve(matches.size());

  for(auto& m : matches) {
    if(m.clippedPoly3D.size() < 3) { // not a snapping matter: left as it is
      kept.push_back(std::move(m));
      continue;
    }

    std::vector<FS_intT> groups;
    std::vector<FSVec3> poly;
    for(const auto& p : m.clippedPoly3D) {
      const FS_intT g = group_[Find(p)];
      if(!groups.empty() && groups.back() == g)
        continue;
      const FS_intT r = rep_[side][g];
      groups.push_back(g);
      poly.push_back(r >= 0 ? points_[r] : p);
    }
    while(groups.size() > 1 && groups.front() == groups.back()) {
      groups.pop_back();
      poly.pop_back();
    }

    // A piece left with less than 3 distinct vertices is the micro-piece itself
    if(std::unordered_set<FS_intT>(groups.begin(), groups.end()).size() < 3) {
      ++dropped;
      continue;
    }
    m.clippedPoly3D = std::move(poly);
    kept.push_back(std::move(m));
  }

  matches = std::move(kept);
  return dropped;
}

_FS_END_NAMESPACE
