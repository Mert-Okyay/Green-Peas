/* This file is auto-generated from gridlink_impl.h.src */
#ifndef DOUBLE_PREC
#define DOUBLE_PREC
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

#include "cellarray_double.h"
#include <inttypes.h>

  extern int get_binsize_double(const double xmin,const double xmax,
                                const double rmax,
                                const int refine_factor,
                                const int max_ncells,
                                double *xbinsize,
                                int *nlattice,
                                const struct config_options *options)  __attribute__((warn_unused_result));

  extern void get_max_min_double(const int64_t ND1, const double * restrict X1, const double * restrict Y1, const double * restrict Z1,
                                 double *min_x, double *min_y, double *min_z, double *max_x, double *max_y, double *max_z);
  

  extern cellarray_double * gridlink_double(const int64_t np,
                                            const double *x,const double *y,const double *z,
                                            const double xmin, const double xmax,
                                            const double ymin, const double ymax,
                                            const double zmin, const double zmax,
                                            const double max_x_size,
                                            const double max_y_size,
                                            const double max_z_size,
                                            const int xbin_refine_factor,
                                            const int ybin_refine_factor,
                                            const int zbin_refine_factor,
                                            int *nlattice_x,
                                            int *nlattice_y,
                                            int *nlattice_z,
                                            const struct config_options *options) __attribute__((warn_unused_result));
  extern void free_cellarray_double(cellarray_double *lattice, const int64_t totncells);
    

  extern cellarray_index_particles_double * gridlink_index_particles_double(const int64_t np,
                                                                            const double *x, const double *y, const double *z, const weight_struct *weights,
                                                                            const double xmin, const double xmax,
                                                                            const double ymin, const double ymax,
                                                                            const double zmin, const double zmax,
                                                                            const double max_x_size,
                                                                            const double max_y_size,
                                                                            const double max_z_size,
                                                                            const int xbin_refine_factor,
                                                                            const int ybin_refine_factor,
                                                                            const int zbin_refine_factor,
                                                                            int *nlattice_x,
                                                                            int *nlattice_y,
                                                                            int *nlattice_z,
                                                                            const struct config_options *options) __attribute__((warn_unused_result));
  extern int assign_ngb_cells_index_particles_double(struct cellarray_index_particles_double *lattice1, struct cellarray_index_particles_double *lattice2,
                                                      const int64_t totncells,
                                                      const int xbin_refine_factor, const int ybin_refine_factor, const int zbin_refine_factor,
                                                      const int nmesh_x, const int nmesh_y, const int nmesh_z,
                                                      const double xdiff, const double ydiff, const double zdiff, 
                                                      const int double_count, const int periodic);
  extern void free_cellarray_index_particles_double(cellarray_index_particles_double *lattice, const int64_t totncells);
  
#ifdef __cplusplus
}
#endif
