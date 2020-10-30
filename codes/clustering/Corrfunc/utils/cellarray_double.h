/* This file is auto-generated from cellarray.h.src */
#ifndef DOUBLE_PREC
#define DOUBLE_PREC
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

#include "weight_defs_double.h"

typedef struct{
  double *x;
  double *y;
  double *z;
  int64_t nelements;//Here the xyz positions will be stored in their individual pointers. More amenable to sorting -> used by wp and xi
} cellarray_double;


typedef struct cellarray_index_particles_double cellarray_index_particles_double;
struct cellarray_index_particles_double{
  int64_t nelements;//Here the xyz positions will be stored in their individual pointers. More amenable to sorting -> used by wp and xi
  int64_t num_ngb;
  double *x;
  double *y;
  double *z;
  weight_struct_double weights;
  cellarray_index_particles_double **ngb_cells;
  double *xwrap;
  double *ywrap;
  double *zwrap;
};

  
#ifdef __cplusplus
}
#endif
