/* This file is auto-generated from weight_functions.h.src */
#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif
// # -*- mode: c -*-
#pragma once

#include "defs.h"
#include "weight_defs_double.h"

#ifdef __AVX__
#include "avx_calls.h"
#endif

#ifdef __SSE4_2__
#include "sse_calls.h"
#endif

#include <stdint.h>

typedef union {
#ifdef __AVX__
    AVX_FLOATS a;
#endif
#ifdef __SSE4_2__
    SSE_FLOATS s;
#endif
    double d;
} weight_union_double;

// Info about a particle pair that we will pass to the weight function
typedef struct
{
    weight_union_double weights0[MAX_NUM_WEIGHTS];
    weight_union_double weights1[MAX_NUM_WEIGHTS];
    weight_union_double dx, dy, dz;
    
    // These will only be present for mock catalogs
    weight_union_double parx, pary, parz;
    
    int64_t num_weights;
} pair_struct_double;

typedef double (*weight_func_t_double)(const pair_struct_double*);
#ifdef __AVX__
typedef AVX_FLOATS (*avx_weight_func_t_double)(const pair_struct_double*);
#endif
#ifdef __SSE4_2__
typedef SSE_FLOATS (*sse_weight_func_t_double)(const pair_struct_double*);
#endif

//////////////////////////////////
// Weighting functions
//////////////////////////////////

/*
 * The pair weight is the product of the particle weights
 */
static inline double pair_product_double(const pair_struct_double *pair){
    return pair->weights0[0].d*pair->weights1[0].d;
}

#ifdef __AVX__
static inline AVX_FLOATS avx_pair_product_double(const pair_struct_double *pair){
    return AVX_MULTIPLY_FLOATS(pair->weights0[0].a, pair->weights1[0].a);
}
#endif

#ifdef __SSE4_2__
static inline SSE_FLOATS sse_pair_product_double(const pair_struct_double *pair){
    return SSE_MULTIPLY_FLOATS(pair->weights0[0].s, pair->weights1[0].s);
}
#endif

//////////////////////////////////
// Utility functions
//////////////////////////////////


/* Gives a pointer to the weight function for the given weighting method
 * and instruction set.
 */
static inline weight_func_t_double get_weight_func_by_method_double(const weight_method_t method){
    switch(method){
        case PAIR_PRODUCT:
            return &pair_product_double;
        default:
        case NONE:
            return NULL;
    }
}

#ifdef __AVX__
static inline avx_weight_func_t_double get_avx_weight_func_by_method_double(const weight_method_t method){
    switch(method){
        case PAIR_PRODUCT:
            return &avx_pair_product_double;
        default:
        case NONE:
            return NULL;
    }
}
#endif

#ifdef __SSE4_2__
static inline sse_weight_func_t_double get_sse_weight_func_by_method_double(const weight_method_t method){
    switch(method){
        case PAIR_PRODUCT:
            return &sse_pair_product_double;
        default:
        case NONE:
            return NULL;
    }
}
#endif
