/* This file is auto-generated from gridlink_impl.c.src */
#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif
// # -*- mode: c -*-
/* File: gridlink_impl.c.src */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "defs.h"
#include "function_precision.h"
#include "utils.h"

#include "gridlink_impl_double.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#define MEMORY_INCREASE_FAC   1.2

int get_binsize_double(const double xmin,const double xmax, const double rmax, const int refine_factor, const int max_ncells, double *xbinsize, int *nlattice, const struct config_options *options)
{
    const double xdiff = (options->periodic && options->boxsize > 0) ? options->boxsize:(xmax-xmin);
    int nmesh=(int)(refine_factor*xdiff/rmax) ;
    nmesh = nmesh < 1 ? 1:nmesh;
    if(options->periodic == 1) {
        if (nmesh<(2*refine_factor+1))  {
            fprintf(stderr,"%s> ERROR:  nlattice = %d is so small that with periodic wrapping the same cells will be counted twice ....exiting\n",__FILE__,nmesh) ;
            fprintf(stderr,"%s> Please reduce Rmax = %"REAL_FORMAT" to be a smaller fraction of the particle distribution region = %"REAL_FORMAT"\n",__FILE__,rmax, xdiff);
            return EXIT_FAILURE;
        }
    }

    if (nmesh>max_ncells)  nmesh=max_ncells;
    *xbinsize = xdiff/nmesh;
    *nlattice = nmesh;

    return EXIT_SUCCESS;
}

void free_cellarray_double(cellarray_double *lattice, const int64_t totncells)
{
    for(int64_t i=0;i<totncells;i++) {
        free(lattice[i].x);
        free(lattice[i].y);
        free(lattice[i].z);
    }
    free(lattice);
}


void free_cellarray_index_particles_double(cellarray_index_particles_double *lattice, const int64_t totncells)
{
    
    for(int64_t i=0;i<totncells;i++){

        free(lattice[i].x);
        free(lattice[i].y);
        free(lattice[i].z);
        for(int w = 0; w < lattice[i].weights.num_weights; w++){
            free(lattice[i].weights.weights[w]);
        }


        /* Might be NULL but free(NULL) is fine*/
        free(lattice[i].xwrap);
        free(lattice[i].ywrap);
        free(lattice[i].zwrap);

        /* Might be NULL but free(NULL) is fine*/
        free(lattice[i].ngb_cells);
    }
    free(lattice);
}


void get_max_min_double(const int64_t ND1, const double * restrict X1, const double * restrict Y1, const double * restrict Z1,
                        double *min_x, double *min_y, double *min_z, double *max_x, double *max_y, double *max_z)
{
    double xmin = *min_x, ymin = *min_y, zmin=*min_z;
    double xmax = *max_x, ymax = *max_y, zmax=*max_z;

    for(int64_t i=0;i<ND1;i++) {
        if(X1[i] < xmin) xmin=X1[i];
        if(Y1[i] < ymin) ymin=Y1[i];
        if(Z1[i] < zmin) zmin=Z1[i];


        if(X1[i] > xmax) xmax=X1[i];
        if(Y1[i] > ymax) ymax=Y1[i];
        if(Z1[i] > zmax) zmax=Z1[i];
    }
    *min_x=xmin;*min_y=ymin;*min_z=zmin;
    *max_x=xmax;*max_y=ymax;*max_z=zmax;
}


cellarray_double * gridlink_double(const int64_t np,
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
                                   const struct config_options *options)
{
    int nmesh_x=0,nmesh_y=0,nmesh_z=0;

    /* Input validation */
    XRETURN(max_x_size > 0.0, NULL, "Minimum separation in X = %"REAL_FORMAT" must be > 0.0\n", max_x_size);
    XRETURN(max_y_size > 0.0, NULL, "Minimum separation in Y = %"REAL_FORMAT" must be > 0.0\n", max_y_size);
    XRETURN(max_z_size > 0.0, NULL, "Minimum separation in Z = %"REAL_FORMAT" must be > 0.0\n", max_z_size);
    XRETURN((xmax-xmin) > 0.0, NULL, "All of the points can not be at the same X position. (xmax-xmin) = %"REAL_FORMAT" must be non-zero\n", (xmax-xmin));
    XRETURN((ymax-ymin) > 0.0, NULL, "All of the points can not be at the same Y position. (ymax-ymin) = %"REAL_FORMAT" must be non-zero\n", (ymax-ymin));
    XRETURN((zmax-zmin) > 0.0, NULL, "All of the points can not be at the same Z position. (zmax-zmin) = %"REAL_FORMAT" must be non-zero\n", (zmax-zmin));
    XRETURN(np > 0, NULL, "Number of points =%"PRId64" must be >0\n", np);
    XRETURN(x != NULL, NULL, "X must be a valid array \n");
    XRETURN(y != NULL, NULL, "Y must be a valid array \n");
    XRETURN(z != NULL, NULL, "Z must be a valid array \n");
    XRETURN(xbin_refine_factor >= 1, NULL, "xbin refine factor = %d must be at least 1\n", xbin_refine_factor);
    XRETURN(ybin_refine_factor >= 1, NULL, "ybin refine factor = %d must be at least 1\n", ybin_refine_factor);
    XRETURN(zbin_refine_factor >= 1, NULL, "zbin refine factor = %d must be at least 1\n", zbin_refine_factor);
    XRETURN(options != NULL, NULL, "Structure containing code options must be a valid address\n");
    
    /* Protect against accidental changes to compile-time (macro) constant */
    XRETURN(MEMORY_INCREASE_FAC >= 1.0, NULL, "Memory increase factor = %lf must be >=1 \n",MEMORY_INCREASE_FAC);


    struct timeval t0;
    if(options->verbose) {
      gettimeofday(&t0,NULL);
    }
    double xbinsize=ZERO, ybinsize=ZERO, zbinsize=ZERO;
    const int xstatus = get_binsize_double(xmin,xmax,max_x_size,xbin_refine_factor, options->max_cells_per_dim, &xbinsize, &nmesh_x, options);
    const int ystatus = get_binsize_double(ymin,ymax,max_y_size,ybin_refine_factor, options->max_cells_per_dim, &ybinsize, &nmesh_y, options);
    const int zstatus = get_binsize_double(zmin,zmax,max_z_size,zbin_refine_factor, options->max_cells_per_dim, &zbinsize, &nmesh_z, options);
    if(xstatus != EXIT_SUCCESS || ystatus != EXIT_SUCCESS || zstatus != EXIT_SUCCESS) {
      return NULL;
    }

    const int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;

    const double xdiff = (options->periodic && options->boxsize > 0) ? options->boxsize:xmax-xmin;
    const double ydiff = (options->periodic && options->boxsize > 0) ? options->boxsize:ymax-ymin;
    const double zdiff = (options->periodic && options->boxsize > 0) ? options->boxsize:zmax-zmin;

    const double cell_volume=xbinsize*ybinsize*zbinsize;
    const double box_volume=xdiff*ydiff*zdiff;
    int64_t expected_n=(int64_t)(np*cell_volume/box_volume*MEMORY_INCREASE_FAC);
    expected_n=expected_n < 2 ? 2:expected_n;

    if(options->verbose) {
      fprintf(stderr,"In %s> Running with [nmesh_x, nmesh_y, nmesh_z]  = %d,%d,%d. ",__FUNCTION__,nmesh_x,nmesh_y,nmesh_z);
    }
    
    cellarray_double *lattice = (cellarray_double *) my_malloc(sizeof(*lattice), totncells);
    int64_t *nallocated = (int64_t *)  my_malloc(sizeof(*nallocated), totncells);
    if(lattice == NULL || nallocated == NULL) {
        free(lattice);free(nallocated);
        return NULL;
    }

    /*
      Allocate memory for each of the fields in cellarray. Since we haven't processed the data yet,
      expected_n is a reasonable guess as to the number of points in the cell.
    */
    for (int64_t index=0;index<totncells;index++) {
        const size_t memsize=sizeof(double);
        lattice[index].x = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
        lattice[index].y = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
        lattice[index].z = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
        if(lattice[index].x == NULL || lattice[index].y == NULL || lattice[index].z == NULL)  {
            for(int64_t j=index;j>=0;j--) {
                free(lattice[j].x);free(lattice[j].y);free(lattice[j].z);
            }
            free(nallocated);free(lattice);
            return NULL;
        }

        lattice[index].nelements=0;
        nallocated[index] = expected_n;
    }

    const double xinv=1.0/xbinsize;
    const double yinv=1.0/ybinsize;
    const double zinv=1.0/zbinsize;

    for (int64_t i=0;i<np;i++)  {
        int ix=(int)((x[i]-xmin)*xinv) ;
        int iy=(int)((y[i]-ymin)*yinv) ;
        int iz=(int)((z[i]-zmin)*zinv) ;
        if (ix>nmesh_x-1)  ix--;    /* this shouldn't happen, but . . . */
        if (iy>nmesh_y-1)  iy--;
        if (iz>nmesh_z-1)  iz--;
        XRETURN(x[i] >= xmin && x[i] <= xmax, NULL, "x[%"PRId64"] = %"REAL_FORMAT" must be within [%"REAL_FORMAT",%"REAL_FORMAT"]\n",
                i, x[i], xmin, xmax);
        XRETURN(y[i] >= ymin && y[i] <= ymax, NULL, "y[%"PRId64"] = %"REAL_FORMAT" must be within [%"REAL_FORMAT",%"REAL_FORMAT"]\n",
                i, y[i], ymin, ymax);
        XRETURN(z[i] >= zmin && z[i] <= zmax, NULL, "z[%"PRId64"] = %"REAL_FORMAT" must be within [%"REAL_FORMAT",%"REAL_FORMAT"]\n",
                i, z[i], zmin, zmax);
        XRETURN(ix >= 0 && ix < nmesh_x, NULL, "ix=%d must be within [0,%d)\n", ix, nmesh_x);
        XRETURN(iy >= 0 && iy < nmesh_y, NULL, "iy=%d must be within [0,%d)\n", iy, nmesh_y);
        XRETURN(iz >= 0 && iz < nmesh_z, NULL, "iz=%d must be within [0,%d)\n", iz, nmesh_z);

        int64_t index = ix*nmesh_y*nmesh_z + iy*nmesh_z + iz;

        if(lattice[index].nelements == nallocated[index]) {
            expected_n = nallocated[index]*MEMORY_INCREASE_FAC;

            //In case expected_n is 1 or MEMORY_INCREASE_FAC is 1.
            //This way, we only increase by a very few particles
            // at a time. Smaller memory footprint
            while(expected_n <= nallocated[index])
                expected_n++;

            const size_t memsize=sizeof(double);
            double *posx=NULL, *posy=NULL, *posz=NULL;
            do {
                posx = my_realloc(lattice[index].x,memsize,expected_n,"lattice.x");
                posy = my_realloc(lattice[index].y,memsize,expected_n,"lattice.y");
                posz = my_realloc(lattice[index].z,memsize,expected_n,"lattice.z");
                if(posx == NULL || posy == NULL || posz == NULL) {
                    expected_n--;
                }
                lattice[index].x = (posx == NULL) ? lattice[index].x:posx;
                lattice[index].y = (posy == NULL) ? lattice[index].y:posy;
                lattice[index].z = (posz == NULL) ? lattice[index].z:posz;
            } while(expected_n > nallocated[index] && (posx == NULL ||
                                                       posy == NULL ||
                                                       posz == NULL));
            
            if(expected_n == nallocated[index]) {
                /*realloc failed. free memory and return */
                fprintf(stderr,"In %s> Reallocation failed,  randomly subsampling the input particle set (currently at %"PRId64" particles) might help\n",
                        __FUNCTION__, np);
                free_cellarray_double(lattice, totncells);
                free(nallocated);
                return NULL;
            }
            nallocated[index] = expected_n;
        }
        XRETURN(lattice[index].nelements < nallocated[index], NULL, 
                ANSI_COLOR_RED"BUG: lattice[%"PRId64"].nelements = %"PRId64" must be less than allocated memory = %"PRId64 ANSI_COLOR_RESET"\n",
                index, lattice[index].nelements, nallocated[index]);

        const int64_t ipos = lattice[index].nelements;
        lattice[index].x[ipos] = x[i];
        lattice[index].y[ipos] = y[i];
        lattice[index].z[ipos] = z[i];
        lattice[index].nelements++;
    }
    free(nallocated);

    //You can free the extra memory reserved by the mallocs by looping over totncells and doing a realloc(lattice[index].x,sizeof(double),lattice[index].nelements,"lattice.x")

    *nlattice_x=nmesh_x;
    *nlattice_y=nmesh_y;
    *nlattice_z=nmesh_z;

    if(options->verbose) {
      struct timeval t1;
      gettimeofday(&t1,NULL);
      fprintf(stderr," Time taken = %7.3lf sec\n",ADD_DIFF_TIME(t0,t1));
    }

    return lattice;
}



/* Need SGLIB to simultaneously sort the particles */
#include "sglib.h"

cellarray_index_particles_double * gridlink_index_particles_double(const int64_t np,
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
                                                                   const struct config_options *options)
{

    int nmesh_x=0,nmesh_y=0,nmesh_z=0;
    
    struct timeval t0;
    if(options->verbose) {
      gettimeofday(&t0,NULL);
    }

    double xbinsize=ZERO, ybinsize=ZERO, zbinsize=ZERO;
    const int xstatus = get_binsize_double(xmin,xmax,max_x_size,xbin_refine_factor, options->max_cells_per_dim, &xbinsize, &nmesh_x, options);
    const int ystatus = get_binsize_double(ymin,ymax,max_y_size,ybin_refine_factor, options->max_cells_per_dim, &ybinsize, &nmesh_y, options);
    const int zstatus = get_binsize_double(zmin,zmax,max_z_size,zbin_refine_factor, options->max_cells_per_dim, &zbinsize, &nmesh_z, options);
    if(xstatus != EXIT_SUCCESS || ystatus != EXIT_SUCCESS || zstatus != EXIT_SUCCESS) {
      fprintf(stderr,"Received xstatus = %d ystatus = %d zstatus = %d. Error\n", xstatus, ystatus, zstatus);
      return NULL;
    }

    const int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;

    const double xdiff = (options->periodic && options->boxsize > 0) ? options->boxsize:xmax-xmin;
    const double ydiff = (options->periodic && options->boxsize > 0) ? options->boxsize:ymax-ymin;
    const double zdiff = (options->periodic && options->boxsize > 0) ? options->boxsize:zmax-zmin;

    const double cell_volume=xbinsize*ybinsize*zbinsize;
    const double box_volume=xdiff*ydiff*zdiff;
    int64_t expected_n=(int64_t)(np*cell_volume/box_volume*MEMORY_INCREASE_FAC);
    expected_n=expected_n < 2 ? 2:expected_n;

    if(options->verbose) {
      fprintf(stderr,"In %s> Running with [nmesh_x, nmesh_y, nmesh_z]  = %d,%d,%d. ",__FUNCTION__,nmesh_x,nmesh_y,nmesh_z);
    }

    cellarray_index_particles_double *lattice  = (cellarray_index_particles_double *) my_malloc(sizeof(*lattice), totncells);
    int64_t *nallocated = (int64_t *)       my_malloc(sizeof(*nallocated)      , totncells);
    if(lattice == NULL || nallocated == NULL) {
        free(lattice);free(nallocated);
        return NULL;
    }
    
    for (int64_t index=0;index<totncells;index++) {
        const size_t memsize=sizeof(double);
        lattice[index].x = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
        lattice[index].y = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
        lattice[index].z = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
        
        // Now do the same for the weights
        lattice[index].weights.num_weights = (weights == NULL) ? 0 : weights->num_weights;
        int w_alloc_status = EXIT_SUCCESS;
        for(int w = 0; w < lattice[index].weights.num_weights; w++){
          lattice[index].weights.weights[w] = (double *) my_malloc(memsize, expected_n);
          if(lattice[index].weights.weights[w] == NULL){
            w_alloc_status = EXIT_FAILURE;
          }
        }
        
        if(lattice[index].x == NULL || lattice[index].y == NULL || lattice[index].z == NULL || w_alloc_status == EXIT_FAILURE) {
            for(int64_t j=index;j>=0;j--) {
                free(lattice[j].x);free(lattice[j].y);free(lattice[j].z);
                for(int w = 0; w < lattice[index].weights.num_weights; w++){
                  free(lattice[index].weights.weights[w]);
                }
            }
            free(nallocated);free(lattice);
            return NULL;
        }

        lattice[index].nelements=0;
        lattice[index].ngb_cells = NULL;
        lattice[index].xwrap = NULL;
        lattice[index].ywrap = NULL;
        lattice[index].zwrap = NULL;
        nallocated[index] = expected_n;
    }

    const double xinv=1.0/xbinsize;
    const double yinv=1.0/ybinsize;
    const double zinv=1.0/zbinsize;

    for (int64_t i=0;i<np;i++)  {
        int ix=(int)((x[i]-xmin)*xinv) ;
        int iy=(int)((y[i]-ymin)*yinv) ;
        int iz=(int)((z[i]-zmin)*zinv) ;

        if (ix>nmesh_x-1)  ix--;    /* this shouldn't happen, but . . . */
        if (iy>nmesh_y-1)  iy--;
        if (iz>nmesh_z-1)  iz--;
        XRETURN(x[i] >= xmin && x[i] <= xmax, NULL, 
               "x[%"PRId64"] = %"REAL_FORMAT" must be within [%"REAL_FORMAT",%"REAL_FORMAT"]\n",
               i, x[i], xmin, xmax);
        XRETURN(y[i] >= ymin && y[i] <= ymax, NULL, 
               "y[%"PRId64"] = %"REAL_FORMAT" must be within [%"REAL_FORMAT",%"REAL_FORMAT"]\n",
               i, y[i], ymin, ymax);
        XRETURN(z[i] >= zmin && z[i] <= zmax, NULL, 
               "z[%"PRId64"] = %"REAL_FORMAT" must be within [%"REAL_FORMAT",%"REAL_FORMAT"]\n",
               i, z[i], zmin, zmax);
        
        XRETURN(ix >= 0 && ix < nmesh_x, NULL, "ix=%d must be within [0,%d)\n", ix, nmesh_x);
        XRETURN(iy >= 0 && iy < nmesh_y, NULL, "iy=%d must be within [0,%d)\n", iy, nmesh_y);
        XRETURN(iz >= 0 && iz < nmesh_z, NULL, "iz=%d must be within [0,%d)\n", iz, nmesh_z);

        const int64_t index = ix*nmesh_y*nmesh_z + iy*nmesh_z + iz;

        if(lattice[index].nelements == nallocated[index]) {
            expected_n = nallocated[index]*MEMORY_INCREASE_FAC;

            //In case expected_n is 1 or MEMORY_INCREASE_FAC is 1.
            //This way, we only increase by a very few particles
            // at a time. Smaller memory footprint
            while(expected_n <= nallocated[index])
                expected_n++;

            const size_t memsize=sizeof(double);
            double *posx=NULL, *posy=NULL, *posz=NULL;
            int w_alloc_status;
            do {
                posx = my_realloc(lattice[index].x,memsize,expected_n,"lattice.x");
                posy = my_realloc(lattice[index].y,memsize,expected_n,"lattice.y");
                posz = my_realloc(lattice[index].z,memsize,expected_n,"lattice.z");
              
                if(posx != NULL){
                  lattice[index].x = posx;
                }
                if(posy != NULL){
                  lattice[index].y = posy;
                }
                if(posz != NULL){
                  lattice[index].z = posz;
                }
                
                w_alloc_status = EXIT_SUCCESS;
                for(int w = 0; w < lattice[index].weights.num_weights; w++){
                  double *newweights = (double *) my_realloc(lattice[index].weights.weights[w], memsize, expected_n, "lattice.weights");
                  if(newweights == NULL){
                    w_alloc_status = EXIT_FAILURE;
                  } else {
                    lattice[index].weights.weights[w] = newweights;
                  }
                }
                
                if(posx == NULL || posy == NULL || posz == NULL || w_alloc_status == EXIT_FAILURE) {
                    expected_n--;
                }
            } while(expected_n > nallocated[index] && (posx == NULL ||
                                                       posy == NULL ||
                                                       posz == NULL ||
                                                       w_alloc_status == EXIT_FAILURE));
            
            if(expected_n == nallocated[index]) {
                /*realloc failed. free memory and return */
                fprintf(stderr,"In %s> Reallocation failed,  randomly subsampling the input particle set (currently at %"PRId64" particles) might help\n",
                        __FUNCTION__, np);
                free_cellarray_index_particles_double(lattice, totncells);
                free(nallocated);
                return NULL;
            }
            nallocated[index] = expected_n;
        }
        XRETURN(lattice[index].nelements < nallocated[index], NULL, 
                ANSI_COLOR_RED"BUG: lattice[%"PRId64"].nelements = %"PRId64" must be less than allocated memory = %"PRId64 ANSI_COLOR_RESET"\n",
                index, lattice[index].nelements, nallocated[index]);

        const int64_t ipos = lattice[index].nelements;
        lattice[index].x[ipos] = x[i];
        lattice[index].y[ipos] = y[i];
        lattice[index].z[ipos] = z[i];
        for(int w = 0; w < lattice[index].weights.num_weights; w++){
            lattice[index].weights.weights[w][ipos] = ((double *)weights->weights[w])[i];
        }
        lattice[index].nelements++;
    }
    free(nallocated);

    /* Do we need to sort the particles in Z ? */
    if(options->sort_on_z) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
        for(int64_t icell=0;icell<totncells;icell++) {
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,X,i,j);        \
                                               SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,Y,i,j);        \
                                               SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,Z,i,j);        \
                                               for(int w = 0; w < first->weights.num_weights; w++){ \
                                                 SGLIB_ARRAY_ELEMENTS_EXCHANGER(double,first->weights.weights[w],i,j);\
                                               }\
                                             }

            const cellarray_index_particles_double *first=&(lattice[icell]);
            if(first->nelements == 0) continue; 
            
            double *X = first->x;
            double *Y = first->y;
            double *Z = first->z;
            
            SGLIB_ARRAY_QUICK_SORT(double, Z, first->nelements, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER
        }
    }


    
    //You can free the extra memory reserved by the mallocs by looping over totncells
    //and doing a realloc(lattice[index].x,sizeof(double),lattice[index].nelements,"lattice.x")
    
    *nlattice_x=nmesh_x;
    *nlattice_y=nmesh_y;
    *nlattice_z=nmesh_z;
    if(options->verbose) {
      struct timeval t1;
      gettimeofday(&t1,NULL);
      fprintf(stderr," Time taken = %7.3lf sec\n",ADD_DIFF_TIME(t0,t1));
    }

    return lattice;
}



int assign_ngb_cells_index_particles_double(struct cellarray_index_particles_double *lattice1, struct cellarray_index_particles_double *lattice2, const int64_t totncells,
                                             const int xbin_refine_factor, const int ybin_refine_factor, const int zbin_refine_factor,
                                             const int nmesh_x, const int nmesh_y, const int nmesh_z,
                                             const double xdiff, const double ydiff, const double zdiff, 
                                             const int autocorr, const int periodic)
{
  const int64_t nx_ngb = 2*xbin_refine_factor + 1;
  const int64_t ny_ngb = 2*ybin_refine_factor + 1;
  const int64_t nz_ngb = 2*zbin_refine_factor + 1;
  const int64_t max_ngb_cells = nx_ngb * ny_ngb * nz_ngb;


  for(int64_t icell=0;icell<totncells;icell++) {
    struct cellarray_index_particles_double *first = &(lattice1[icell]);
    if(first->nelements == 0) continue;
    const int iz = icell % nmesh_z;
    const int ix = icell / (nmesh_y * nmesh_z );
    const int iy = (icell - iz - ix*nmesh_z*nmesh_y)/nmesh_z;
    XRETURN(icell == (ix * nmesh_y * nmesh_z + iy * nmesh_z + (int64_t) iz), EXIT_FAILURE,
            ANSI_COLOR_RED"BUG: Index reconstruction is wrong. icell = %"PRId64" reconstructed index = %"PRId64 ANSI_COLOR_RESET"\n",
            icell, (ix * nmesh_y * nmesh_z + iy * nmesh_z + (int64_t) iz));
    
    first->num_ngb = 0;
    if(periodic == 1) {
      first->xwrap = my_malloc(sizeof(*(first->xwrap)), max_ngb_cells);
      first->ywrap = my_malloc(sizeof(*(first->ywrap)), max_ngb_cells);
      first->zwrap = my_malloc(sizeof(*(first->zwrap)), max_ngb_cells);
      if(first->xwrap == NULL || first->ywrap == NULL || first->zwrap == NULL)  {
          return EXIT_FAILURE;
      }
    } else {
        first->xwrap = NULL;
        first->ywrap = NULL;
        first->zwrap = NULL;
    }
    first->ngb_cells = my_malloc(sizeof(*(first->ngb_cells)) , max_ngb_cells);
    if(first->ngb_cells == NULL) {
        return EXIT_FAILURE;
    }
    
    for(int iix=-xbin_refine_factor;iix<=xbin_refine_factor;iix++){
      const int periodic_ix = (ix + iix + nmesh_x) % nmesh_x;
      const int non_periodic_ix = ix + iix;
      const int iiix = (periodic == 1) ? periodic_ix:non_periodic_ix;
      if(iiix < 0 || iiix >= nmesh_x) continue;
      const double off_xwrap = ((ix + iix) >= 0) && ((ix + iix) < nmesh_x) ? 0.0: ((ix+iix) < 0 ? xdiff:-xdiff);
      
      for(int iiy=-ybin_refine_factor;iiy<=ybin_refine_factor;iiy++) {
        const int periodic_iy = (iy + iiy + nmesh_y) % nmesh_y;
        const int non_periodic_iy = iy + iiy;
        const int iiiy = (periodic == 1) ? periodic_iy:non_periodic_iy;
        if(iiiy < 0 || iiiy >= nmesh_y) continue;
        const double off_ywrap = ((iy + iiy) >= 0) && ((iy + iiy) < nmesh_y) ? 0.0: ((iy+iiy) < 0 ? ydiff:-ydiff);
          const int start_iz = -zbin_refine_factor;
          for(int64_t iiz=start_iz;iiz<=zbin_refine_factor;iiz++){
            const int periodic_iz = (iz + iiz + nmesh_z) % nmesh_z;
            const int non_periodic_iz = iz + iiz;
            const int iiiz = (periodic == 1) ? periodic_iz:non_periodic_iz;
            if(iiiz < 0 || iiiz >= nmesh_z) continue;
            
            const double off_zwrap = ((iz + iiz) >= 0) && ((iz + iiz) < nmesh_z) ? 0.0: ((iz+iiz) < 0 ? zdiff:-zdiff);
            const int64_t icell2 = iiiz + (int64_t) nmesh_z*iiiy + nmesh_z*nmesh_y*iiix;
            
            //For cases where we are not double-counting (i.e., wp and xi), the same-cell
            //must always be evaluated. In all other cases, (i.e., where double-counting is occurring)
            //is used, include that in the ngb_cells! The interface is a lot cleaner in the double-counting
            //kernels in that case! Also, if the second cell has no particles, then just skip it
            if((autocorr == 1 && icell2 >= icell) || lattice2[icell2].nelements == 0) {
                continue;
            }
            const int64_t ngb_index = first->num_ngb;
            XRETURN(ngb_index < max_ngb_cells, EXIT_FAILURE,
                    "ngb index = %"PRId64" should be less than max_ngb = %"PRId64"\n", ngb_index, max_ngb_cells);
            first->ngb_cells[ngb_index] = &(lattice2[icell2]);
            
            //Note the xwrap/ywraps do not have memory allocated for them in the
            //non-periodic case. 
            if(periodic == 1) {
              first->xwrap[ngb_index] = off_xwrap;
              first->ywrap[ngb_index] = off_ywrap;
              first->zwrap[ngb_index] = off_zwrap;
            } 
            first->num_ngb++;
          }
      }
    }
  }
  
  return EXIT_SUCCESS;
}    

