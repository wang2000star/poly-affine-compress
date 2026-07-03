#ifndef WALSH_H
#define WALSH_H

#include "types.h"

std::vector<WalshInfo> compute_walsh_correlations(const TruthTable& tt, int n_threads);

#endif
