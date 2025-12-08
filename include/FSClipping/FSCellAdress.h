#ifndef FSCELLADRESS_H
#define FSCELLADRESS_H

#include "FSClac.h"
#include "FSCommon.h"
#include "FSConfig.h"
#include "FSTypes.h"

struct FSCellAddress {
  /// The process where the cell resides.
  FS_intT procID;
  /// The cellID identifying the cell (which is only valid on procID).
  FS_intT cellID;

  /// Copy constructor.
  /**
   * @param other The object to copy construct from.
   */
  FSCellAddress(const FSCellAddress &other)
      : procID(other.procID), cellID(other.cellID) {}

  /// Set from data constructor.
  /**
   * @param procID The process index.
   * @param cellID The cell index.
   */
  FSCellAddress(const FS_intT procID, const FS_intT cellID)
      : procID(procID), cellID(cellID) {}

  /// Copy-assignment operator from another cell address.
  /**
   * @param rhs The cell address whose value we want to assign to `this`.
   */
  void operator=(const FSCellAddress &rhs) {
    procID = rhs.procID;
    cellID = rhs.cellID;
  }
  /// Equality comparison operator for cell address.
  /**
   * @param rhs The right-hand side argument of the comparison-operator.
   * @return true if `this` and rhs are equal, false otherwise.
   */
  bool operator==(const FSCellAddress &rhs) const {
    return std::tie(procID, cellID) == std::tie(rhs.procID, rhs.cellID);
  }

  /// Inequality comparison operator for cell address.
  /**
   * @param rhs The right-hand side argument of the comparison-operator.
   * @return true if `this` and rhs are different, false otherwise.
   */
  bool operator!=(const FSCellAddress &rhs) const {
    return std::tie(procID, cellID) != std::tie(rhs.procID, rhs.cellID);
  }

  /// Less-than comparison operator for cell address.
  /**
   * The less-than comparison operator defines an ordering first by domain, then
   * cell-index on that domain.
   *
   * @param rhs The right-hand side argument of the comparison-operator.
   * @return true if `this` is less than rhs; else false.
   */
  bool operator<(const FSCellAddress &rhs) const {
    return std::tie(procID, cellID) < std::tie(rhs.procID, rhs.cellID);
  }

  /// Get the size of the object when buffered for FSDM MPI communication.
  /**
   * @param clac The FSDM MPI communications handler object.
   * @return the size of the object when buffered for FSDM MPI communication.
   */
  typename FSClac::sizeT GetBufSize(FSClac &clac) const {
    return clac.GetBufSize<FS_intT>() + clac.GetBufSize<FS_intT>();
  }

  /// Pack object into FSClac active MPI communication buffer.
  /**
   * @param clac The FSDM MPI communications handler object.
   */
  void Pack(FSClac &clac) {
    clac.Pack(&procID);
    clac.Pack(&cellID);
  }

  /// Unpack object from FSClac active MPI communication buffer.
  /**
   * @param clac The FSDM MPI communications handler object.
   */
  void Unpack(FSClac &clac) {
    clac.Unpack(&procID);
    clac.Unpack(&cellID);
  }
};

#endif