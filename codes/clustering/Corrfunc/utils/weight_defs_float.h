/* This file is auto-generated from weight_defs.h.src */
#ifdef DOUBLE_PREC
#undef DOUBLE_PREC
#endif
// # -*- mode: c -*-
#pragma once

#include <stdint.h>

typedef struct
{
    float *weights[MAX_NUM_WEIGHTS];  // This will be of shape weights[num_weights][num_particles]
    int64_t num_weights;
} weight_struct_float;
