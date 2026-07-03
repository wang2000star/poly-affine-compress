#ifndef ANF_H
#define ANF_H

#include "types.h"

int64_t count_T(const uint64_t* data, int n);
RawANFInfo compute_raw_anf_info(TruthTable& tt);

#endif
