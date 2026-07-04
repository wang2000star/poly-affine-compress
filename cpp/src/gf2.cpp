#include "gf2.h"
#include <cstring>

int gf2_rank(const uint32_t* rows, int m, int n) {
    uint32_t tmp[32];
    for (int i = 0; i < m; i++) tmp[i] = rows[i];
    int rank = 0;
    for (int col = n - 1; col >= 0 && rank < m; col--) {
        int pivot = -1;
        for (int i = rank; i < m; i++) {
            if ((tmp[i] >> col) & 1) { pivot = i; break; }
        }
        if (pivot < 0) continue;
        uint32_t t = tmp[rank]; tmp[rank] = tmp[pivot]; tmp[pivot] = t;
        for (int i = 0; i < m; i++) {
            if (i != rank && ((tmp[i] >> col) & 1))
                tmp[i] ^= tmp[rank];
        }
        rank++;
    }
    return rank;
}

bool gf2_invert(const uint32_t* M, int n, uint32_t* M_inv) {
    uint32_t tmp[32][32];
    uint32_t inv[32][32];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tmp[i][j] = (M[i] >> (n - 1 - j)) & 1;
            inv[i][j] = (i == j) ? 1 : 0;
        }
    }
    for (int col = 0; col < n; col++) {
        int pivot = -1;
        for (int i = col; i < n; i++) {
            if (tmp[i][col]) { pivot = i; break; }
        }
        if (pivot < 0) return false;
        for (int j = 0; j < n; j++) {
            uint32_t t = tmp[col][j]; tmp[col][j] = tmp[pivot][j]; tmp[pivot][j] = t;
            t = inv[col][j]; inv[col][j] = inv[pivot][j]; inv[pivot][j] = t;
        }
        for (int i = 0; i < n; i++) {
            if (i != col && tmp[i][col]) {
                for (int j = 0; j < n; j++) {
                    tmp[i][j] ^= tmp[col][j];
                    inv[i][j] ^= inv[col][j];
                }
            }
        }
    }
    for (int i = 0; i < n; i++) {
        M_inv[i] = 0;
        for (int j = 0; j < n; j++)
            M_inv[i] |= (inv[i][j] << (n - 1 - j));
    }
    return true;
}

bool gf2_right_inverse(const uint32_t* M, int m, int n, uint32_t* P) {
    uint32_t M_work[32], A[32];
    for (int i = 0; i < m; i++) {
        M_work[i] = M[i];
        A[i] = (1u << i);
    }

    int pivot_cols[32];
    int rank = 0;

    for (int col = 0; col < n && rank < m; col++) {
        int pivot = -1;
        for (int i = rank; i < m; i++) {
            if ((M_work[i] >> col) & 1) { pivot = i; break; }
        }
        if (pivot < 0) continue;

        uint32_t t = M_work[rank]; M_work[rank] = M_work[pivot]; M_work[pivot] = t;
        t = A[rank]; A[rank] = A[pivot]; A[pivot] = t;
        pivot_cols[rank] = col;

        for (int i = 0; i < m; i++) {
            if (i != rank && ((M_work[i] >> col) & 1)) {
                M_work[i] ^= M_work[rank];
                A[i] ^= A[rank];
            }
        }
        rank++;
    }

    if (rank < m) return false;

    memset(P, 0, n * sizeof(uint32_t));
    for (int i = 0; i < m; i++)
        P[pivot_cols[i]] = A[i];
    return true;
}

int gf2_kernel_basis(const uint32_t* M, int m, int n, int* pivot_cols,
                     uint32_t* kernel, uint32_t* M_work) {
    uint32_t work[32];
    if (M_work) {
        for (int i = 0; i < m; i++) work[i] = M_work[i];
    } else {
        for (int i = 0; i < m; i++) work[i] = M[i];
    }

    int rank = 0;
    int pivots[32];

    for (int col = 0; col < n && rank < m; col++) {
        int pivot = -1;
        for (int i = rank; i < m; i++) {
            if ((work[i] >> col) & 1) { pivot = i; break; }
        }
        if (pivot < 0) continue;

        uint32_t t = work[rank]; work[rank] = work[pivot]; work[pivot] = t;
        pivots[rank] = col;

        for (int i = 0; i < m; i++) {
            if (i != rank && ((work[i] >> col) & 1))
                work[i] ^= work[rank];
        }
        rank++;
    }

    if (pivot_cols) {
        for (int i = 0; i < rank; i++) pivot_cols[i] = pivots[i];
    }

    // For each free column, construct kernel basis vector
    int k_idx = 0;
    int pivot_idx = 0;
    for (int f = 0; f < n; f++) {
        if (pivot_idx < rank && pivots[pivot_idx] == f) {
            pivot_idx++;
            continue;
        }
        // Free column f: construct kernel vector
        uint32_t v = (1u << f);  // bit f set
        for (int r = 0; r < rank; r++) {
            if ((work[r] >> f) & 1)
                v |= (1u << pivots[r]);
        }
        if (kernel) kernel[k_idx] = v;
        k_idx++;
    }

    return k_idx;
}
