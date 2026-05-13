#ifndef FSFACEEXCHANGE_H
#define FSFACEEXCHANGE_H

#include "FSClipping/FSClippingFace.h"
#include <FSClac.h>
#include <vector>

_FS_BEGIN_NAMESPACE

// Point-to-point exchange of FSClippingFace geometry via the FSClac buffer
// protocol. Only vertex geometry is transmitted; FSDM topology (_faceFSDM)
// is left null on the receiving side, which is acceptable since
// ComputeMatches guards all _faceFSDM accesses.
namespace FSFaceExchange {

// Pack face geometry into clac's send buffer and transmit to destProc.
// destProc must call Receive(clac, thisProc).
void Send(FSClac& clac, FS_intT destProc, const std::vector<FSClippingFace>& faces);

// Receive and unpack face geometry transmitted by sourceProc via Send.
std::vector<FSClippingFace> Receive(FSClac& clac, FS_intT sourceProc);

} // namespace FSFaceExchange

_FS_END_NAMESPACE

#endif // FSFACEEXCHANGE_H
