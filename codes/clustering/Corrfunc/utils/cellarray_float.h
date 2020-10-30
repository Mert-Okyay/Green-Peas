/* This file is auto-generated from cellarray.h.src */
#ifdef DOUBLE_PREC
#undef DOUBLE_PREC
#endif
// # -*- mode: c -*-
/* File: cellarray.h.src */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdint.h>

#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "weight_defs_float.h"

typedef struct{
  float *x;
  float *y;
  float *z;
  int64_t nelements;//Here the xyz positions will be stored in their individual pointers. More amenable to sorting -> used by wp and xi
} cellarray_float;


typedef struct cellarray_index_particles_float cellarray_index_particles_float;
struct cellarray_index_particles_float{
  int64_t nelements;//Here the xyz positions will be stored in their individual pointers. More amenable to sorting -> used by wp and xi
  int64_t num_ngb;
  float *x;
  float *y;
  float *z;
  weight_struct_float weights;
  cellarray_index_particles_float **ngb_cells;
  float *xwrap;
  float *ywrap;
  float *zwrap;
};

  
#ifdef __cplusplus
}
#endif
