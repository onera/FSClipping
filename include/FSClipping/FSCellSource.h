#ifndef FSCELLSOURCE_H
#define FSCELLSOURCE_H


#include <FSCellPool.h>
#include <FSMesh.h>

_FS_BEGIN_NAMESPACE

class FSCellSource
{
public:
  FSCellSource(FSMeshEnums::CellType t, const FSCellPool& cellPool) : cellType_(t), cellPool_(cellPool){};

  const FSMeshEnums::CellType& Type() const noexcept { return cellType_; }
  bool IsPoly2D() const;
  bool IsPoly3D() const;

  FS_intT NNodes(FS_intT c) const;
  FS_intT Node(FS_intT c, FS_intT n) const;

private:
  FSMeshEnums::CellType cellType_;
  const FSCellPool& cellPool_;
};

std::vector<FSCellSource> BuildCellSources(const FSMesh& mesh);

_FS_END_NAMESPACE

#endif