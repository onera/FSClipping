#include "FSClipping/FSVertexSnap.h"
#include "FSClipping/FSFaceMatcher.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
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

FS_floatT DistanceToOutline(const FSVec3& p, const std::vector<FSVec3>& poly)
{
  FS_floatT best = kInf;
  for(std::size_t i = 0; i < poly.size(); ++i) {
    const FSVec3& a = poly[i];
    const FSVec3& b = poly[(i + 1) % poly.size()];
    const FSVec3 ab = b - a;
    const FS_floatT len2 = ab.InnerProduct(ab);
    const FS_floatT t = len2 > 0.0 ? std::clamp((p - a).InnerProduct(ab) / len2, 0.0, 1.0) : 0.0;
    const FSVec3 q(a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]);
    best = std::min(best, q.Distance(p));
  }
  return best;
}

FS_intT Root(std::vector<FS_intT>& parent, FS_intT i)
{
  while(parent[i] != i) {
    parent[i] = parent[parent[i]];
    i = parent[i];
  }
  return i;
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

static std::array<std::uint64_t, 3> KeyOf(const FSVec3& p)
{
  auto bits = [](FS_floatT x) { return std::bit_cast<std::uint64_t>(x == 0.0 ? 0.0 : x); };
  return {bits(p[0]), bits(p[1]), bits(p[2])};
}

FS_intT FSVertexSnap::AddPoint(const FSVec3& p)
{
  const auto [it, inserted] = indexOf_.try_emplace(KeyOf(p), FS_intT(points_.size()));
  if(inserted) {
    points_.push_back(p);
    for(FS_intT s = 0; s < 2; ++s) {
      isOriginal_[s].push_back(false);
      outlineDist_[s].push_back(kInf);
    }
  }
  return it->second;
}

FS_intT FSVertexSnap::Find(const FSVec3& p) const
{
  const auto it = indexOf_.find(KeyOf(p));
  if(it == indexOf_.end())
    Fail("vertex missing from the raw matches", p);
  return it->second;
}

FSVertexSnap::FSVertexSnap(const std::vector<FSClippingFace>& faces1, const std::vector<FSClippingFace>& faces2,
                           const std::vector<FSFaceMatch>& rawMatches, FS_floatT snapTol)
{
  // 1 - Collect the points: original vertices, then the raw pieces in the frame of each side
  const std::array<const std::vector<FSClippingFace>*, 2> faces{&faces1, &faces2};
  std::array<std::unordered_map<FS_intT, const FSClippingFace*>, 2> faceOf;
  for(FS_intT s = 0; s < 2; ++s)
    for(const auto& f : *faces[s]) {
      faceOf[s].emplace(f.faceIndex(), &f);
      for(const auto& v : f.vertices()) {
        const FS_intT i = AddPoint(v);
        isOriginal_[s][i] = true;
        outlineDist_[s][i] = 0.0;
      }
    }

  auto addPiece = [&](FS_intT s, FS_intT face, const std::vector<FSVec3>& poly) {
    const auto it = faceOf[s].find(face);
    for(const auto& p : poly) {
      const FS_intT i = AddPoint(p);
      if(it != faceOf[s].end())
        outlineDist_[s][i] = std::min(outlineDist_[s][i], DistanceToOutline(p, it->second->vertices()));
    }
  };
  for(const auto& m : rawMatches) {
    addPiece(0, m.face1, m.clippedPoly3D);
    addPiece(1, m.face2, m.clippedPoly3D_face2);
  }

  // 2 - Groups: transitive closure of the pairs closer than snapTol
  const FS_intT n = FS_intT(points_.size());
  std::vector<FS_intT> parent(n);
  std::iota(parent.begin(), parent.end(), 0);
  std::unordered_map<CellKey, std::vector<FS_intT>, CellKeyHash> grid;
  for(FS_intT i = 0; i < n; ++i) {
    const CellKey c = CellOf(points_[i], snapTol);
    for(std::int64_t dx = -1; dx <= 1; ++dx)
      for(std::int64_t dy = -1; dy <= 1; ++dy)
        for(std::int64_t dz = -1; dz <= 1; ++dz) {
          const auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
          if(it == grid.end())
            continue;
          for(const FS_intT j : it->second)
            if(points_[i].Distance(points_[j]) <= snapTol)
              parent[Root(parent, i)] = Root(parent, j);
        }
    grid[c].push_back(i);
  }

  group_.assign(n, -1);
  std::unordered_map<FS_intT, FS_intT> groupOfRoot;
  for(FS_intT i = 0; i < n; ++i)
    group_[i] = groupOfRoot.try_emplace(Root(parent, i), FS_intT(groupOfRoot.size())).first->second;

  // 3 - One representative per side: its original node, else the point nearest to its face outline
  const FS_intT nGroups = FS_intT(groupOfRoot.size());
  for(FS_intT s = 0; s < 2; ++s) {
    rep_[s].assign(nGroups, -1);
    std::vector<bool> repIsOriginal(nGroups, false);
    for(FS_intT i = 0; i < n; ++i) {
      const FS_intT g = group_[i];
      FS_intT& r = rep_[s][g];
      if(isOriginal_[s][i]) {
        // Merging two nodes of one mesh would collapse one of its edges
        if(repIsOriginal[g])
          Fail("two nodes of the same mesh closer than the snap tolerance", points_[i]);
        r = i;
        repIsOriginal[g] = true;
      } else if(!repIsOriginal[g] && outlineDist_[s][i] < kInf && (r < 0 || outlineDist_[s][i] < outlineDist_[s][r]))
        r = i;
    }
  }
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
