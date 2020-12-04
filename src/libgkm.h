/* File: libgkm.h
 *
 * Description: header file for gkm module for libsvm
 *
 * Copyright (C) 2015 Dongwon Lee
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef LIBSVM_GKM_H_INCLUDED
#define LIBSVM_GKM_H_INCLUDED

//#include "libsvm.h"
#define LOGGER_ID 0
#define LOGGER_FORMAT "%l %d %t: %m\n"

#define MAX_ALPHABET_SIZE 4 /* base ACGT, DON'T CHANGE THIS PARAMETER! */
#define MAX_ALPHABET_SIZE_SQ 16 // MAX_ALPHABET_SIZE*MAX_ALPHABET_SIZE
#define MAX_SEQ_LENGTH 2048
#define MMCNT_LOOKUPTAB_WIDTH 8

#ifdef __cplusplus 
extern "C" {
#endif

typedef struct _KmerTree KmerTree;
typedef struct _KmerTreeLeaf KmerTreeLeaf;
typedef struct _NodeMismatchCount NodeMismatchCount;
typedef struct _KmerTreeCoef KmerTreeCoef;
typedef struct _KmerTreeLeafData KmerTreeLeafData;
typedef struct _gkm_parameter gkm_parameter;
typedef struct _gkm_data gkm_data;
typedef struct _svm_problem svm_problem;

enum { GKM, EST_FULL, EST_TRUNC, EST_TRUNC_RBF, EST_TRUNC_PW, EST_TRUNC_PW_RBF}; /* kernel_type */

struct _gkm_parameter {
    int kernel_type;
    int L;
    int k;
    int d;
    uint8_t M;
    double H;
    double gamma;
};

struct _gkm_data {
    char *sid;
    int seqid;
    int label;
    int seqlen;
    uint8_t *seq;
    uint8_t *seq_rc;
    uint8_t *wt;
    uint8_t *wt_rc;
    int *kmerids;
    int *kmerids_rc;
    char *seq_string;

    double sqnorm;
};

struct _svm_problem
{
    int l;
    double *y;
    gkm_data **x;
};

struct _KmerTreeLeafData {
    int seqid;
    int wt;
};

/* simple stack */
struct _KmerTreeLeaf {
    int count;
    int capacity;
    KmerTreeLeafData *data;
};

struct _KmerTree {
    int depth; //the same as kmer-length
    int node_count; //internal node only
    int leaf_count; //leaf node only
    int *node;
    KmerTreeLeaf *leaf;
};

struct _KmerTreeCoef {
    int depth; //the same as kmer-length
    int node_count; //internal node only
    int leaf_count; //leaf node only
    double *coef_sum;
};

void gkmkernel_init(gkm_parameter *param);
void gkmkernel_destroy();
void gkmkernel_set_num_threads(int n);

gkm_data* gkmkernel_new_object(char *seq, char *sid, int seqid);
void gkmkernel_delete_object(gkm_data* d);
void gkmkernel_free_object(gkm_data* d);

double gkmkernel_kernelfunc(const gkm_data *da, const gkm_data *db);
double* gkmkernel_kernelfunc_batch(const gkm_data *da, const gkm_data **db_array, const int n, double *res);

void gkmkernel_init_problems(gkm_data **x, int n);
void gkmkernel_destroy_problems();
void gkmkernel_swap_index(int i, int j);
void gkmkernel_update_index();
double* gkmkernel_kernelfunc_batch_all(const int a, const int start, const int end, double *res);

void gkmkernel_destroy_sv();
double* gkmkernel_kernelfunc_batch_sv(const gkm_data *d, double *res);


#ifdef __cplusplus
}
#endif

#endif
