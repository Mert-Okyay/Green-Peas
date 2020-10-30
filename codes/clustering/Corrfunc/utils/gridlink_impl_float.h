/* This file is auto-generated from gridlink_impl.h.src */
#ifdef DOUBLE_PREC
#undef DOUBLE_PREC
#endif
// # -*- mode: c -*-
/* File: gridlink_impl.h.src */
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

#include "cellarray_float.h"
#include <inttypes.h>

  extern int get_binsize_float(const float xmin,const float xmax,
                                const float rmax,
                                const int refine_factor,
                                const int max_ncells,
                                float *xbinsize,
                                int *nlattice,
                                const struct config_options *options)  __attribute__((warn_unused_result));

  extern void get_max_min_float(const int64_t ND1, const float * restrict X1, const float * restrict Y1, const float * restrict Z1,
                                 float *min_x, float *min_y, float *min_z, float *max_x, float *max_y, float *max_z);
  

  extern cellarray_float * gridlink_float(const int64_t np,
                                            const float *x,const float *y,const float *z,
                                            const float xmin, const float xmax,
                                            const float ymin, const float ymax,
                                            const float zmin, const float zmax,
                                            const float max_x_size,
                                            const float max_y_size,
                                            const float max_z_size,
                                            const int xbin_refine_factor,
                                            const int ybin_refine_factor,
                                            const int zbin_refine_factor,
                                            int *nlattice_x,
                                            int *nlattice_y,
                                            int *nlattice_z,
                                            const struct config_options *options) __attribute__((warn_unused_result));
  extern void free_cellarray_float(cellarray_float *lattice, const int64_t totncells);
    

  extern cellarray_index_particles_float * gridlink_index_particles_float(const int64_t np,
                                                                            const float *x, const float *y, const float *z, const weight_struct *weights,
                                                                            const float xmin, const float xmax,
                                                                            const float ymin, const float ymax,
                                                                            const float zmin, const float zmax,
                                                                            const float max_x_size,
                                                                            const float max_y_size,
                                                                            const float max_z_size,
                                                                            const int xbin_refine_factor,
                                                                            const int ybin_refine_factor,
                                                                            const int zbin_refine_factor,
                                                                            int *nlattice_x,
                                                                            int *nlattice_y,
                                                                            int *nlattice_z,
                                                                            const struct config_options *options) __attribute__((warn_unused_result));
  extern int assign_ngb_cells_index_particles_float(struct cellarray_index_particles_float *lattice1, struct cellarray_index_particles_float *lattice2,
                                                      const int64_t totncells,
                                                      const int xbin_refine_factor, const int ybin_refine_factor, const int zbin_refine_factor,
                                                      const int nmesh_x, const int nmesh_y, const int nmesh_z,
                                                      const float xdiff, const float ydiff, const float zdiff, 
                                                      const int double_count, const int periodic);
  extern void free_cellarray_index_particles_float(cellarray_index_particles_float *lattice, const int64_t totncells);
  
#ifdef __cplusplus
}
#endif
