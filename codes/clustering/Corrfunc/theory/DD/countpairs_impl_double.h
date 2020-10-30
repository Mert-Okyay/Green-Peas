/* This file is auto-generated from countpairs_impl.h.src */
#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif
// # -*- mode: c -*-
/* File: countpairs_impl.h.src */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#include "defs.h"
#include "weight_defs_double.h"
#include <inttypes.h>

#include "countpairs.h"  /* For definition of results_countpairs */

    extern void interrupt_handler_countpairs_double(int signo);
    
    typedef int (*countpairs_func_ptr_double)(const int64_t N0, double *x0, double *y0, double *z0, const weight_struct_double *weights0,
                                             const int64_t N1, double *x1, double *y1, double *z1, const weight_struct_double *weights1,
                                             const int same_cell,
                                             const double sqr_rpmax, const double sqr_rpmin, const int nbin, const double *rupp_sqr, const double rpmax,
                                             const double off_xwrap, const double off_ywrap, const double off_zwrap,
                                             double *src_rpavg, uint64_t *src_npairs,
                                             double *src_weightavg, const weight_method_t weight_method);
  
    extern countpairs_func_ptr_double countpairs_driver_double(const struct config_options *options) __attribute__((warn_unused_result));
    
    extern int countpairs_double(const int64_t ND1, double *X1, double *Y1, double  *Z1,
                                 const int64_t ND2, double *X2, double *Y2, double  *Z2,
                                 const int numthreads,
                                 const int autocorr,
                                 const char *binfile,
                                 results_countpairs *results,
                                 struct config_options *options,
                                 struct extra_options *extra);

#ifdef __cplusplus
}
#endif
