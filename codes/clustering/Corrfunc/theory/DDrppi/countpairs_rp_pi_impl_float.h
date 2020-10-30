/* This file is auto-generated from countpairs_rp_pi_impl.h.src */
#ifdef DOUBLE_PREC
#undef DOUBLE_PREC
#endif
// # -*- mode: c -*-
/* File: countpairs_rp_pi_impl.h.src */
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

#include "defs.h" //for struct config_options 
#include "weight_defs_float.h"
#include <inttypes.h> //for uint64_t

#include "countpairs_rp_pi.h"//for struct results_countpairs_rp_pi

    extern void interrupt_handler_countpairs_rp_pi_float(int signo);
    
    typedef int (*countpairs_rp_pi_func_ptr_float)(const int64_t N0, float *x0, float *y0, float *z0, const weight_struct_float *weights0,
                                                    const int64_t N1, float *x1, float *y1, float *z1, const weight_struct_float *weights1, const int same_cell,
                                                    const float sqr_rpmax, const float sqr_rpmin, const int nbin, const int npibin,
                                                    const float *rupp_sqr, const float pimax,
                                                    const float off_xwrap, const float off_ywrap, const float off_zwrap,
                                                    float *src_rpavg, uint64_t *src_npairs,
                                                    float *src_weightavg, const weight_method_t weight_method);

    
    extern countpairs_rp_pi_func_ptr_float countpairs_rp_pi_driver_float(const struct config_options *options) __attribute__((warn_unused_result));

    extern int countpairs_rp_pi_float(const int64_t ND1, float *X1, float *Y1, float *Z1,
                                       const int64_t ND2, float *X2, float *Y2, float *Z2,
                                       const int numthreads,
                                       const int autocorr,
                                       const char *binfile,
                                       const float pimax,
                                       results_countpairs_rp_pi *results,
                                       struct config_options *options,
                                       struct extra_options *extra);

#ifdef __cplusplus
}
#endif
