#include "FSClipping/FSCellSource.h"
#include <FSMeshData.h>

_FS_BEGIN_NAMESPACE

bool FSCellSource::IsPoly2D() const
{
  bool isPoly2D = false;
  if(FSMeshEnums::IsUnstructSurfaceCellType(cellType_) && FSMeshEnums::IsPolyCellType(cellType_))
    isPoly2D = true;
  return isPoly2D;
};

bool FSCellSource::IsPoly3D() const
{
  bool isPoly3D = false;
  if(FSMeshEnums::IsUnstructVolumeCellType(cellType_) && FSMeshEnums::IsPolyCellType(cellType_))
    isPoly3D = true;
  return isPoly3D;
};

FS_intT FSCellSource::NNodes(FS_intT c) const { return FSCellInfo::NNodes(cellType_, cellPool_, c); }

FS_intT FSCellSource::Node(FS_intT c, FS_intT n) const { return FSCellInfo::GetCellNode(cellType_, cellPool_, c, n); }

std::vector<FSCellSource> BuildCellSources(FSMesh& mesh)
{
  const auto& cellTypes = mesh.GetCellTypes();
  std::vector<FSCellSource> cellSource;

  for(const auto& t : mesh.GetCellTypes()) {
    const auto& cellPool = mesh.GetMeshData()->GetUnstructCells().GetCellPool(t);
    cellSource.emplace_back(t, *cellPool);
  }

  return cellSource;
}
_FS_END_NAMESPACE