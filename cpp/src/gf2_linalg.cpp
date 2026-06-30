#include "gf2_linalg.h"
#include <algorithm>
#include <cassert>
#include <iostream>

int gf2_rank(const GF2Matrix& M, int n) {
    GF2Matrix rows = M;
    int rank = 0;
    for (int col = 0; col < n && rank < (int)rows.size(); ++col) {
        // Find pivot row
        int pivot = -1;
        for (int r = rank; r < (int)rows.size(); ++r) {
            if ((rows[r] >> col) & 1) {
                pivot = r;
                break;
            }
        }
        if (pivot < 0) continue;

        std::swap(rows[rank], rows[pivot]);

        // Eliminate other rows
        uint64_t mask = rows[rank];
        for (int r = 0; r < (int)rows.size(); ++r) {
            if (r != rank && ((rows[r] >> col) & 1)) {
                rows[r] ^= mask;
            }
        }
        ++rank;
    }
    return rank;
}

GF2Matrix gf2_invert(const GF2Matrix& M, int n) {
    // Augment [M | I]
    GF2Matrix aug = M;
    GF2Matrix inv(n, 0);
    for (int i = 0; i < n; ++i) {
        inv[i] = (uint64_t)1 << i;
    }

    for (int col = 0; col < n; ++col) {
        // Find pivot
        int pivot = -1;
        for (int r = col; r < n; ++r) {
            if ((aug[r] >> col) & 1) {
                pivot = r;
                break;
            }
        }
        if (pivot < 0) return {}; // not invertible

        std::swap(aug[col], aug[pivot]);
        std::swap(inv[col], inv[pivot]);

        uint64_t mask = aug[col];
        uint64_t inv_mask = inv[col];

        for (int r = 0; r < n; ++r) {
            if (r != col && ((aug[r] >> col) & 1)) {
                aug[r] ^= mask;
                inv[r] ^= inv_mask;
            }
        }
    }
    return inv;
}

GF2Matrix gf2_extend_to_invertible(const GF2Matrix& M, int m, int n) {
    assert(m < n);
    GF2Matrix rows = M;
    rows.resize(n, 0);

    // Gaussian elimination on the first m rows to find pivot columns
    std::vector<bool> pivot_col(n, false);
    int rank = 0;
    for (int col = 0; col < n && rank < m; ++col) {
        int pivot = -1;
        for (int r = rank; r < m; ++r) {
            if ((rows[r] >> col) & 1) {
                pivot = r;
                break;
            }
        }
        if (pivot < 0) continue;
        std::swap(rows[rank], rows[pivot]);
        pivot_col[col] = true;
        uint64_t mask = rows[rank];
        for (int r = 0; r < m; ++r) {
            if (r != rank && ((rows[r] >> col) & 1)) {
                rows[r] ^= mask;
            }
        }
        ++rank;
    }

    // Fill missing rows with unit vectors for non-pivot columns
    int slot = m;
    for (int col = 0; col < n; ++col) {
        if (!pivot_col[col]) {
            rows[slot++] = (uint64_t)1 << col;
        }
    }

    return rows;
}

GF2Matrix gf2_extend_columns_to_invertible(const GF2Matrix& M, int m, int n) {
    assert(m > n);
    // Gaussian elimination on M to find pivot rows
    GF2Matrix tmp = M;
    std::vector<int> perm(m);
    for (int i = 0; i < m; ++i) perm[i] = i;

    int rank = 0;
    for (int col = 0; col < n && rank < m; ++col) {
        int pivot = -1;
        for (int r = rank; r < m; ++r) {
            if ((tmp[r] >> col) & 1) {
                pivot = r;
                break;
            }
        }
        if (pivot < 0) continue;
        std::swap(tmp[rank], tmp[pivot]);
        std::swap(perm[rank], perm[pivot]);
        uint64_t mask = tmp[rank];
        for (int r = 0; r < m; ++r) {
            if (r != rank && ((tmp[r] >> col) & 1)) {
                tmp[r] ^= mask;
            }
        }
        ++rank;
    }
    // rank should equal n (full column rank). Non-pivot rows at perm[n..m-1].

    // Build m×m extended matrix by adding unit vector columns for non-pivot rows
    GF2Matrix result(m, 0);
    for (int r = 0; r < m; ++r) {
        result[r] = M[r];  // first n bits from original M
    }
    for (int j = 0; j < m - n; ++j) {
        int np_row = perm[n + j];
        // Column (n+j) has a 1 at position np_row (and 0 elsewhere)
        result[np_row] |= (uint64_t)1 << (n + j);
    }

    return result;
}

uint64_t gf2_mat_vec_mul(const GF2Matrix& M, uint64_t v) {
    uint64_t result = 0;
    for (int i = 0; i < (int)M.size(); ++i) {
        if (__builtin_parityll(M[i] & v)) {
            result |= (uint64_t)1 << i;
        }
    }
    return result;
}

void gf2_print_matrix(const GF2Matrix& M, int rows, int cols) {
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            std::cout << ((M[r] >> c) & 1);
        }
        std::cout << "\n";
    }
}
