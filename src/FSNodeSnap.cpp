#include "FSClipping/FSNodeSnap.h"

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <unordered_map>

_FS_BEGIN_NAMESPACE

namespace {

// Bitwise identity of a point: the face vertices are exact copies of the mesh coordinates
using ExactKey = std::array<std::uint64_t, 3>;

ExactKey KeyOf(const FSVec3& p)
{
  auto bits = [](FS_floatT x) { return std::bit_cast<std::uint64_t>(x == 0.0 ? 0.0 : x); };
  return {bits(p[0]), bits(p[1]), bits(p[2])};
}

struct ExactKeyHash {
  std::size_t operator()(const ExactKey& k) const noexcept
  {
    std::size_t seed = 0;
    for(const auto v : k)
      seed ^= std::hash<std::uint64_t>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

// Cell of a grid of step snapTol: two points closer than snapTol are in neighbouring cells
using CellKey = std::array<std::int64_t, 3>;

CellKey CellOf(const FSVec3& p, FS_floatT step)
{
  return {static_cast<std::int64_t>(std::floor(p[0] / step)), static_cast<std::int64_t>(std::floor(p[1] / step)),
          static_cast<std::int64_t>(std::floor(p[2] / step))};
}

struct CellKeyHash {
  std::size_t operator()(const CellKey& k) const noexcept
  {
    std::size_t seed = 0;
    for(const auto v : k)
      seed ^= std::hash<std::int64_t>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

std::vector<FSVec3> DistinctVertices(const std::vector<FSClippingFace>& faces)
{
  std::unordered_map<ExactKey, FS_intT, ExactKeyHash> seen;
  std::vector<FSVec3> points;
  for(const auto& f : faces)
    for(const auto& v : f.vertices())
      if(seen.try_emplace(KeyOf(v), FS_intT(points.size())).second)
        points.push_back(v);
  return points;
}

void Fail(const char* what, const FSVec3& p)
{
  FSString msg = "FSNodeSnap: ";
  msg += what;
  msg += " near (";
  msg.Add(p[0]);
  msg += ", ";
  msg.Add(p[1]);
  msg += ", ";
  msg.Add(p[2]);
  msg += "), decrease the snap tolerance";
  FSError.SetAndPrintAndExit(msg.c_str());
}

} // namespace

std::vector<FSNodeMove> FSNodeSnap::Compute(const std::vector<FSClippingFace>& target,
                                            const std::vector<FSClippingFace>& moving, FS_floatT snapTol)
{
  std::vector<FSNodeMove> moves;
  if(snapTol <= 0.0)
    return moves;

  const std::vector<FSVec3> targets = DistinctVertices(target);
  std::unordered_map<CellKey, std::vector<FS_intT>, CellKeyHash> grid;
  for(FS_intT t = 0; t < FS_intT(targets.size()); ++t)
    grid[CellOf(targets[t], snapTol)].push_back(t);

  std::vector<bool> taken(targets.size(), false);
  for(const FSVec3& p : DistinctVertices(moving)) {
    const CellKey c = CellOf(p, snapTol);
    FS_intT found = -1;
    for(std::int64_t dx = -1; dx <= 1; ++dx)
      for(std::int64_t dy = -1; dy <= 1; ++dy)
        for(std::int64_t dz = -1; dz <= 1; ++dz) {
          const auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
          if(it == grid.end())
            continue;
          for(const FS_intT t : it->second) {
            if(targets[t].Distance(p) > snapTol)
              continue;
            // Two candidates: snapping would collapse an edge of the target mesh
            if(found >= 0)
              Fail("several target nodes within the snap tolerance", p);
            found = t;
          }
        }

    if(found < 0 || KeyOf(targets[found]) == KeyOf(p))
      continue;
    if(taken[found])
      Fail("two moving nodes snapped onto the same target node", targets[found]);
    taken[found] = true;
    moves.push_back({p, targets[found]});
  }
  return moves;
}

void FSNodeSnap::Apply(const std::vector<FSNodeMove>& moves, std::vector<FSClippingFace>& faces)
{
  if(moves.empty())
    return;
  std::unordered_map<ExactKey, FSVec3, ExactKeyHash> to;
  for(const auto& m : moves)
    to.emplace(KeyOf(m.from), m.to);

  for(auto& f : faces) {
    std::vector<FSVec3> verts = f.vertices();
    bool moved = false;
    for(auto& v : verts) {
      const auto it = to.find(KeyOf(v));
      if(it == to.end())
        continue;
      v = it->second;
      moved = true;
    }
    if(moved)
      f.setVertices(verts);
  }
}

FS_intT FSNodeSnap::Apply(const std::vector<FSNodeMove>& moves, std::vector<FSVec3>& points)
{
  if(moves.empty())
    return 0;
  std::unordered_map<ExactKey, FSVec3, ExactKeyHash> to;
  for(const auto& m : moves)
    to.emplace(KeyOf(m.from), m.to);

  FS_intT nMoved = 0;
  for(auto& p : points) {
    const auto it = to.find(KeyOf(p));
    if(it == to.end())
      continue;
    p = it->second;
    ++nMoved;
  }
  return nMoved;
}

_FS_END_NAMESPACE
