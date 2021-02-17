/*
The eXtended Keccak Code Package (XKCP)
https://github.com/XKCP/XKCP

Keccak, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van
Assche.

Implementation by the designers, hereby denoted as "the implementer".

For more information, feedback or questions, please refer to the Keccak Team
website: https://keccak.team/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "KeccakSponge.h"

#include "KeccakP-1600-SnP.h"

#define prefix KeccakWidth1600
#define SnP KeccakP1600
#define SnP_width 1600
#define SnP_Permute KeccakP1600_Permute_24rounds
#if defined(KeccakF1600_FastLoop_supported)
#define SnP_FastLoop_Absorb KeccakF1600_FastLoop_Absorb
#endif
#include "KeccakSponge.inc"
#undef prefix
#undef SnP
#undef SnP_width
#undef SnP_Permute
#undef SnP_FastLoop_Absorb

#define prefix KeccakWidth1600_12rounds
#define SnP KeccakP1600
#define SnP_width 1600
#define SnP_Permute KeccakP1600_Permute_12rounds
#if defined(KeccakP1600_12rounds_FastLoop_supported)
#define SnP_FastLoop_Absorb KeccakP1600_12rounds_FastLoop_Absorb
#endif
#include "KeccakSponge.inc"
#undef prefix
#undef SnP
#undef SnP_width
#undef SnP_Permute
#undef SnP_FastLoop_Absorb
