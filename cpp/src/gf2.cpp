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
