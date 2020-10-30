/* This file is auto-generated from weight_functions.h.src */
#ifdef DOUBLE_PREC
#undef DOUBLE_PREC
#endif
// # -*- mode: c -*-
#pragma once

#include "defs.h"
#include "weight_defs_float.h"

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
    float d;
} weight_union_float;

// Info about a particle pair that we will pass to the weight function
typedef struct
{
    weight_union_float weights0[MAX_NUM_WEIGHTS];
    weight_union_float weights1[MAX_NUM_WEIGHTS];
    weight_union_float dx, dy, dz;
    
    // These will only be present for mock catalogs
    weight_union_float parx, pary, parz;
    
    int64_t num_weights;
} pair_struct_float;

typedef float (*weight_func_t_float)(const pair_struct_float*);
#ifdef __AVX__
typedef AVX_FLOATS (*avx_weight_func_t_float)(const pair_struct_float*);
#endif
#ifdef __SSE4_2__
typedef SSE_FLOATS (*sse_weight_func_t_float)(const pair_struct_float*);
#endif

//////////////////////////////////
// Weighting functions
//////////////////////////////////

/*
 * The pair weight is the product of the particle weights
 */
static inline float pair_product_float(const pair_struct_float *pair){
    return pair->weights0[0].d*pair->weights1[0].d;
}

#ifdef __AVX__
static inline AVX_FLOATS avx_pair_product_float(const pair_struct_float *pair){
    return AVX_MULTIPLY_FLOATS(pair->weights0[0].a, pair->weights1[0].a);
}
#endif

#ifdef __SSE4_2__
static inline SSE_FLOATS sse_pair_product_float(const pair_struct_float *pair){
    return SSE_MULTIPLY_FLOATS(pair->weights0[0].s, pair->weights1[0].s);
}
#endif

//////////////////////////////////
// Utility functions
//////////////////////////////////


/* Gives a pointer to the weight function for the given weighting method
 * and instruction set.
 */
static inline weight_func_t_float get_weight_func_by_method_float(const weight_method_t method){
    switch(method){
        case PAIR_PRODUCT:
            return &pair_product_float;
        default:
        case NONE:
            return NULL;
    }
}

#ifdef __AVX__
static inline avx_weight_func_t_float get_avx_weight_func_by_method_float(const weight_method_t method){
    switch(method){
        case PAIR_PRODUCT:
            return &avx_pair_product_float;
        default:
        case NONE:
            return NULL;
    }
}
#endif

#ifdef __SSE4_2__
static inline sse_weight_func_t_float get_sse_weight_func_by_method_float(const weight_method_t method){
    switch(method){
        case PAIR_PRODUCT:
            return &sse_pair_product_float;
        default:
        case NONE:
            return NULL;
    }
}
#endif
