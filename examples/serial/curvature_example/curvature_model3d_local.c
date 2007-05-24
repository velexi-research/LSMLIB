// System headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* LSMLIB headers */
#include "lsm_level_set_evolution3d_local.h"
#include "lsm_spatial_derivatives3d_local.h"
#include "lsm_utilities3d_local.h"
#include "lsm_tvd_runge_kutta3d_local.h"
#include "lsm_reinitialization3d_local.h"
#include "lsm_geometry3d.h"
#include "lsm_localization3d.h"

/* LSMLIB Serial package headers */
#include "lsm_boundary_conditions.h"
#include "lsm_macros.h"

/* Local headers */
#include "curvature_model_top.h"
#include "curvature_model3d.h"
#include "curvature_model3d_local.h"

#define DT_MIN_TO_CORRECT 1e-5
#define DT_MIN            0.001

/* limits for average gradient of phi within the narrow band: it
   shouldn't be too far away from 1.0 */
#define AVE_GRAD_PHI_MIN 0.9
#define AVE_GRAD_PHI_MAX 1.1

static unsigned char mark_gb=127, mark_D1=126, mark_D2=125, mark_fb=124;
 
 	
/* 
* Localization (narrow banding) implementation follows
* Peng/Merriman/Osher/Zhao/Kang paper "A PDE-Based Fast Local Level Set Method"
* Journal of Computational Physics, 1999. 
*/

void curvatureModelMedium3dLocalMainLoop(
     Options          *options,
     LSM_DataArrays   *data_arrays,
     Grid             *grid,
     FILE             *fp_out)
{
  double   cfl_number = 0.5;
  
  /* time variables */
  double   t, dt, dt_sub, max_H, dt_corr;
  double   tplot, dt_min, dt_max;
  
  double   max_abs_err, eps, eps_stop;
 
  double   zero = 0.0;
  double   vel_n, vol_phi, vol_max, vol_phi_prev, rel_vol_diff;
  int      i, nx, nxy;  
  
  int      bdry_location_idx = 9; /* extrapolate all boundaries */
  
  int      OUTER_STEP, INNER_STEP, TOTAL_STEP;
  int      reinit_steps, last_reinit_step, ave_reinit_steps;
  
  
  /* writing shortcuts */
  Grid             *g = grid;
  LSM_DataArrays   *d = data_arrays;
  Options          *o = options;
   
  /* variables specific for localization */
  double   beta, gamma;
  int      nlo_index, nhi_index, level; 
  
  double   frac_nb, last_reinit_time, grad_phi_ave;  
  int      nb_level0, nb_level1, nb_level2;
  int      reinit_trigger;
  
  int      nlo_index_outer, nhi_index_outer;
  int      n_outer, change_sgn;
  int      n_lo_copy[6], n_hi_copy[6];
  int      change_sgn_steps, grad_phi_ave_steps;
   
  t = 0;
  /* every TPLOT time period we evaluate max. abs. error as well as
  *  reinitalize the function as needed.
  */
  max_abs_err = 1000.0;
  
  /* stopping criterion - somewhat arbitrary,
    modify EMAX_STOP in the accompanying .h file
  */
  eps_stop = EMAX_STOP *(g->dx)[0]; //stopping criterion
  if( options->print_details)
  {    
    fprintf(fp_out,"\nTPLOT %g eps_stop %g set internally\n",TPLOT,eps_stop);
    fprintf(fp_out,"Simulation continues until given time tmax is reached\n"); 
    fprintf(fp_out,"or max.abs.error for phi(:,t) - phi(:,t-TPLOT)");
    fprintf(fp_out," is less than eps_stop.\n");
    fprintf(fp_out,"-----------------------------------------------------\n");
    fprintf(fp_out,"\nEach TPLOT time, we report on:");
    fprintf(fp_out,"\nvol_phi - the volume occupied by the neg. level set phase");
    fprintf(fp_out,"\nvol_frac - the fraction of pore space (vol_max)");
    fprintf(fp_out," occupied by the neg. level set phase");
    fprintf(fp_out,"\nrel_vol_diff - the relative volume difference btw two time steps");   
    fprintf(fp_out,"\n");
  }
  
  /* correction for time spacing due to parabolic (curvature) term */
  dt_corr = 1.0/((g->dx)[0]*(g->dx)[0]) + 1.0/((g->dx)[1]*(g->dx)[1]) + 
            1.0/((g->dx)[2]*(g->dx)[2]);
  dt_corr *= 2.0*o->b;
  
  /* this eps is suggested for Heaviside function in Fedkiw/Osher book */
  eps = 1.5*(g->dx[0]);
  nx = (g->grid_dims_ghostbox)[0];
  nxy = (g->grid_dims_ghostbox)[0]*(g->grid_dims_ghostbox)[1];
    
  /* compute volume of the pore space */
  LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_max,
	        d->mask,
		&(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		&(g->klo_gb), &(g->khi_gb),
		&(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		&(g->klo_fb), &(g->khi_fb),
		&(g->dx[0]),&(g->dx[1]),&(g->dx[2]),
		&eps);
  
  /* compute volume of the fluid */
  LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi,
	        d->phi,
		&(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		&(g->klo_gb), &(g->khi_gb),
		&(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		&(g->klo_fb), &(g->khi_fb),
		&(g->dx[0]),&(g->dx[1]),&(g->dx[2]),
		&eps);	 
  	
  tplot = (TPLOT < o->tmax ) ? TPLOT : o->tmax;
   		
  dt_min = 100.0; dt_max = 0;
  
  /* localization: parameters for the 2nd order ENO or 3rd order WENO 
     defining the narrow band width */
  beta = 2*(g->dx)[0]; gamma = 4*(g->dx)[0];
  /* will need narrow band level 0,1,2,3 */
  level = 3; 
  
  /* localization: reinitialize globally so T0 can be set */
  reinitializeMedium3d(d,g,o,gamma + g->dx[0]); 	 

  /* localization - allocated number of index points array elements */
  nlo_index = 0; 
  nhi_index = g->num_gridpts - 1;
  nlo_index_outer = 0;
  nhi_index_outer = d->num_alloc_index_outer_pts-1;
  
  OUTER_STEP = 0; INNER_STEP = 0; TOTAL_STEP = 0;
  last_reinit_step = 0; last_reinit_time = 0;
  reinit_steps = change_sgn_steps = grad_phi_ave_steps = 0;
  ave_reinit_steps = 0;
  
  while( (t < o->tmax)  && (max_abs_err > eps_stop) && (vol_phi > eps_stop))
  {  /* outer loop - the code is set up to output some error information
     *  volume fractions etc. as the computation progresses
     */
    OUTER_STEP++;
    dt_sub = 0;
    
    COPY_DATA(d->phi_prev,d->phi,g)
    vol_phi_prev = vol_phi; 
    
    while( dt_sub < tplot )
    { /* inner loop */
      INNER_STEP++;
      TOTAL_STEP++;
      
       /* localization : determine T0 */	   
      LSM3D_DETERMINE_NARROW_BAND(d->phi,
           &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
           &(g->klo_gb), &(g->khi_gb),
	   d->narrow_band,
	   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
           &(g->klo_gb), &(g->khi_gb),
	   d->index_x, d->index_y, d->index_z,
	   &nlo_index, &nhi_index,	   
	   d->n_lo,d->n_hi,
	   d->index_outer_pts,
	   &nlo_index_outer, &nhi_index_outer,
	   &(d->nlo_outer_plus),  &(d->nhi_outer_plus),
	   &(d->nlo_outer_minus), &(d->nhi_outer_minus),
           &gamma,&beta,&level);
	   
     
      /* mark boundary layers in narrow_band array 
      *  These layer marks to be used in Fortran functions for checking if the
      *	 point is in the correct fill box.
     */     	   	   
      LSM3D_MARK_NARROW_BAND_BOUNDARY_LAYER(d->narrow_band,
       	   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
           &(g->klo_gb), &(g->khi_gb),
	   &(g->ilo_D2_fb), &(g->ihi_D2_fb), &(g->jlo_D2_fb), &(g->jhi_D2_fb),
           &(g->klo_D2_fb), &(g->khi_D2_fb),
	   &mark_D2);
	   
      LSM3D_MARK_NARROW_BAND_BOUNDARY_LAYER(d->narrow_band,
       	   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
           &(g->klo_gb), &(g->khi_gb),
	   &(g->ilo_D1_fb), &(g->ihi_D1_fb), &(g->jlo_D1_fb), &(g->jhi_D1_fb),
           &(g->klo_D1_fb), &(g->khi_D1_fb),
	   &mark_D1);
	   
      LSM3D_MARK_NARROW_BAND_BOUNDARY_LAYER(d->narrow_band,
       	   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
           &(g->klo_gb), &(g->khi_gb),
	   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
           &(g->klo_gb), &(g->khi_gb),
	   &mark_gb);	   	   
	   
      LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL(d->lse_rhs,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
                    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0]);
     
      if(o->a > 0)
      {  
         /* Compute upwinding gradient approximations */ 
          LSM3D_HJ_ENO2_LOCAL(d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
                    d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->D1,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->D2,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    &(d->n_lo)[1],&(d->n_hi)[1],
		    &(d->n_lo)[2],&(d->n_hi)[2],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb,&mark_D1,&mark_D2); 
	 
	 vel_n = o->a;
	 
	 LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL(d->lse_rhs,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &vel_n,
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],		    
                    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb);
                  
	 /* figure out time spacing for hyperbolic term */
	 LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL(&dt,&vel_n,
		    d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_minus, d->phi_y_minus,  d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),		    
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),&cfl_number,
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],		    
                    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),	
		    &mark_fb);  
      }
      else dt = tplot;
     
      if( o->b > 0)
      {
	/* Compute derivatives needed for curvature */	
	LSM3D_CENTRAL_GRAD_ORDER2_LOCAL(d->phi_x, d->phi_y, d->phi_z,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z, 
		    &(d->n_lo)[0],&(d->n_hi)[1],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_D1);	
	LSM3D_CENTRAL_GRAD_ORDER2_LOCAL(d->phi_xx, d->phi_xy, d->phi_xz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z, 
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_D2);		
	LSM3D_CENTRAL_GRAD_ORDER2_LOCAL(d->phi_xy, d->phi_yy, d->phi_yz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_y,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_D2);	
	LSM3D_CENTRAL_GRAD_ORDER2_LOCAL(d->phi_xz, d->phi_yz, d->phi_zz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_z,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_D2);
		    
	LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL(d->lse_rhs,
	            &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x,d->phi_y,d->phi_z,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_xx,d->phi_xy,d->phi_xz,
		    d->phi_yy,d->phi_yz,d->phi_zz,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(o->b),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0], &(d->n_hi)[0],
                    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb);
		    	
	/* correct dt due to parabolic (curvature) term */
        if( o->a > 0 )
	   max_H = cfl_number / dt;
	else 
	   max_H = 0;
	   
	dt =  cfl_number / (max_H + dt_corr);
      }	    

      if( dt < DT_MIN_TO_CORRECT ) dt = DT_MIN; 
       
      if(dt_sub + dt > tplot) 
      {
	  dt = tplot - dt_sub;
      }      

      /* collect info on max. and min. time spacing */
      if(dt > dt_max) dt_max = dt;
      if(dt < dt_min) dt_min = dt;
      
      /* localization: modify equation by a cut-off function */
      LSM3D_MULTIPLY_CUT_OFF_LSE_RHS_LOCAL(d->phi, d->lse_rhs,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb,
		    &beta,&gamma);
      
      LSM3D_TVD_RK2_STAGE1_LOCAL(d->phi_stage1,
                   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->phi,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->lse_rhs,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &dt,
		   d->index_x, d->index_y, d->index_z,
		   &(d->n_lo)[0],&(d->n_hi)[0],
		   d->narrow_band,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &mark_fb);	

       /* boundary conditions */
       signedLinearExtrapolationBC(d->phi_stage1,g,bdry_location_idx);
         
      /* masking enforced so that the interface stays within pore space */
      if(o->do_mask) IMPOSE_MASK_LOCAL(d->phi_stage1,d->mask,d->phi_stage1,g,d);       

      LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL(d->lse_rhs,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
                    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[2]);  
      
      if(o->a)
      {
	  LSM3D_HJ_ENO2_LOCAL(d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
                    d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_stage1,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->D1,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->D2,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    &(d->n_lo)[1],&(d->n_hi)[1],
		    &(d->n_lo)[2],&(d->n_hi)[2],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb,&mark_D1,&mark_D2);   
		        
	 LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL(d->lse_rhs,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &vel_n,
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],		    
                    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb);	    
      }
      
      if( o->b )
      {       
	LSM3D_CENTRAL_GRAD_ORDER2_LOCAL(d->phi_x,d->phi_y,d->phi_z,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_stage1,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[1],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_D1);		
	LSM3D_CENTRAL_GRAD_ORDER2_LOCAL(d->phi_xx,d->phi_xy,d->phi_xz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_D2);		
	LSM3D_CENTRAL_GRAD_ORDER2_LOCAL(d->phi_xy,d->phi_yy,d->phi_yz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_y,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_D2);	
	LSM3D_CENTRAL_GRAD_ORDER2_LOCAL(d->phi_xz, d->phi_yz, d->phi_zz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_z,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_D2);

	LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL(d->lse_rhs,
	            &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x,d->phi_y,d->phi_z,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_xx,d->phi_xy,d->phi_xz,
		    d->phi_yy,d->phi_yz,d->phi_zz,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(o->b),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0], &(d->n_hi)[0],
                    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb);
      }     
     
      /* localization: modify equation by a cut-off function */
      LSM3D_MULTIPLY_CUT_OFF_LSE_RHS_LOCAL(d->phi_stage1, d->lse_rhs,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],		    
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb,&beta,&gamma);
		    
      LSM3D_TVD_RK2_STAGE2_LOCAL(d->phi_next,
                   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->phi_stage1,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->phi,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->lse_rhs,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &dt,
		   d->index_x, d->index_y, d->index_z,
		   &(d->n_lo)[0],&(d->n_hi)[0],
		   d->narrow_band,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &mark_fb);
        
      /* boundary conditions */
       signedLinearExtrapolationBC(d->phi_next,g,bdry_location_idx);	 
           
      /* masking enforced so that the interface stays within pore space */
      if(o->do_mask) IMPOSE_MASK_LOCAL(d->phi,d->mask,d->phi_next,g,d)
      else           COPY_DATA(d->phi,d->phi_next,g)
           
      /* localization : check if the sign of the level set function
                       changes in the outer layer of the narrow band
		       where  gamma > |phi| >= beta  */
      LSM3D_CHECK_OUTER_NARROW_BAND_LAYER(&change_sgn,
            d->phi,
            &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
            &(g->klo_gb), &(g->khi_gb),
            d->index_x, d->index_y, d->index_z,
            &(d->n_lo)[0],&(d->n_hi)[0],
            d->index_outer_pts,
	    &nlo_index_outer, &nhi_index_outer,
	    &(d->nlo_outer_plus),  &(d->nhi_outer_plus),
	    &(d->nlo_outer_minus), &(d->nhi_outer_minus));      
      
      if(change_sgn)
      {  /* if the sign changed, the interface is close to the narrow
            band border so reinitialization needs to be triggerred */
         reinit_trigger = 1;
	 change_sgn_steps++; 
      }	 
      else
      {
         /* localization : compute the average value of the norm
	 of the gradient */	  
	 
         LSM3D_COMPUTE_AVE_GRAD_PHI_LOCAL(&grad_phi_ave,
	      d->phi,
              &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
              &(g->klo_gb), &(g->khi_gb),	 
              &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
              d->index_x, d->index_y, d->index_z,
              &(d->n_lo)[0],&(d->n_hi)[0],
              d->narrow_band,
              &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	      &(g->klo_gb), &(g->khi_gb),
	      &mark_fb);
	//printf("\n grad_phi_ave %g", grad_phi_ave); fflush(stdout);     
	if(( grad_phi_ave < AVE_GRAD_PHI_MIN ) || 
	   ( grad_phi_ave > AVE_GRAD_PHI_MAX ))
	{
	     /* if the norm "considerably" deviates from 1.0, trigger
	        reinitialization */
	    reinit_trigger = 1;
	    grad_phi_ave_steps++;
	}
	else
	{  /* reinitialize if haven't done so for long time */ 
	   if( t + dt_sub - last_reinit_time >= TPLOT )
	       reinit_trigger = 1;
	   else 
	       reinit_trigger = 0;   
	}            
      }
      
      /* localization: reinitalize if needed */
      if(reinit_trigger)
       {
           /* gather some stats on how often we reinitialize */
	   last_reinit_time = t;
           ave_reinit_steps += TOTAL_STEP - last_reinit_step;
           last_reinit_step = TOTAL_STEP;
	   reinit_steps++;	   
	   
           /* N0 for reinitialization purposes is tube T0 plus its
	     first neighbors; essentially level 0 and level 1 narrow band */
	     
	   for(i = 0; i < 6; i++)
	   {
	     n_lo_copy[i] = d->n_lo[i];    n_hi_copy[i] = d->n_hi[i];
	   }
	   
	   /* shift limits in order to reinitialize on wider narrow band */
	   d->n_hi[0] = d->n_hi[1];
	   d->n_lo[1] = d->n_lo[2];   d->n_hi[1] = d->n_hi[2];
	   d->n_lo[2] = d->n_lo[3];   d->n_hi[2] = d->n_hi[3];
	         
           reinitializeMedium3dLocal(d,g,o,gamma + 2*g->dx[0]);
	   
	   /* copy old limit values back */
	   for(i = 0; i < 6; i++)
	   {
	      d->n_lo[i] = n_lo_copy[i];  d->n_hi[i] = n_hi_copy[i];
	   }
       }
    
      dt_sub = dt_sub + dt;
   } /* inner loop */
  
   t = t + dt_sub;   
  
   /* compute max abs error only in the narrow band */
   LSM3D_MAX_NORM_DIFF_LOCAL(&max_abs_err,d->phi,
            &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	    &(g->klo_gb), &(g->khi_gb),
	    d->phi_prev,
	    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	    &(g->klo_gb), &(g->khi_gb),
	    d->index_x, d->index_y, d->index_z,
            &(d->n_lo)[0],&(d->n_hi)[0],
            d->narrow_band,
	    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
            &(g->klo_gb), &(g->khi_gb),	
            &mark_fb);
   LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi,
	    d->phi,
            &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	    &(g->klo_gb), &(g->khi_gb),
            &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
            &(g->klo_gb), &(g->khi_gb),
            &(g->dx[0]),&(g->dx[1]),&(g->dx[2]),
	    &eps);		    
    
   printf("Time interval [%g,%g], max. abs. error %g\n", t-tplot,t,max_abs_err);
   fprintf(fp_out,"Time interval [%g,%g], max. abs. error %g\n",t-tplot,t,
                                                                  max_abs_err);
   /* relative difference in volume */
   rel_vol_diff = fabs(vol_phi_prev - vol_phi)/vol_phi_prev;		
   fprintf(fp_out," rel_vol_diff %g vol_phi %g vol_frac %g\n",
                                 rel_vol_diff,vol_phi,vol_phi/vol_max);
   fprintf(fp_out," dt_min %g dt_max %g\n",dt_min,dt_max);
   fflush(stdout); fflush(fp_out);
   
   /* checking the number of points on the narrow band */
   nb_level0 = (d->n_hi)[0] - (d->n_lo)[0] + 1;
   nb_level1 = (d->n_hi)[1] - (d->n_lo)[1] + 1;
   nb_level2 = (d->n_hi)[2] - (d->n_lo)[2] + 1;
   frac_nb = (nb_level0 + nb_level1 + nb_level2)/(double)g->num_gridpts;
   fprintf(fp_out,"narrow band level0 %8d all levels %d total frac %g\n",
             nb_level0,nb_level0+nb_level1+nb_level2,frac_nb);
   /* 
   n_outer = d->nhi_outer_plus -  d->nlo_outer_plus + 1 + 
             d->nhi_outer_minus - d->nlo_outer_minus + 1;	     
   fprintf(fp_out," n_outer %d\n",n_outer);
   */	     
   
	     
   fflush(stdout); fflush(fp_out);								  										  
  } /* outer loop */
   
  ave_reinit_steps =ceil( (double)(ave_reinit_steps) / (double)(reinit_steps) );
  fprintf(fp_out,"\nTotal steps %d   Reinit. steps %d  (change sign %d, grad_phi_ave %d)",
         TOTAL_STEP, reinit_steps,change_sgn_steps,grad_phi_ave_steps);
  fprintf(fp_out,"\nReinitialized on average every %d steps.\n",ave_reinit_steps); 
  
}


/* 
*  reinitializeMedium3dLocal() reinitializes the level set function using the 
*  second order accuracy ENO and TVD RK routines.
*  The computation is performed locally (within the narrow band).
*  
*  Arguments:
*   data_arrays  - LSMLIB Serial package data arrays structure
*   grid         - LSMLIB Serial package Grid structure
*   options      - (local) Options structure; the only element used is 'do_mask'
*   tmax_r -     Maximal running time for reinitialization. Note that the normal
*               velocity in the reinitialization level set equation is 1. This 
*	        time is hence equal to the distance from the interface within 
*	        which the level set function will be replaced by a signed 
*	        distance function.     
*/
   
void reinitializeMedium3dLocal(
     LSM_DataArrays *data_arrays,
     Grid           *grid,
     Options        *options,
     double         tmax_r)
{   
    double cfl_number = 0.5;
    double t_r, dt_r;
    
    int    use_phi0_for_sign = 0;
    int    nx, nxy;
    int    bdry_location_idx = 9; /* all boundaries */
  
     /* writing shortcuts */
    Grid             *g = grid;
    LSM_DataArrays   *d = data_arrays;
    Options          *o = options;
    
    nx =  (g->grid_dims_ghostbox)[0];
    nxy = (g->grid_dims_ghostbox)[0]*(g->grid_dims_ghostbox)[1];           
    
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
      
    COPY_DATA(d->phi0,d->phi,g)
    
    while(t_r < tmax_r )
    {
       LSM3D_HJ_ENO2_LOCAL(d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb),&(g->khi_gb),
                    d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb),&(g->khi_gb),
		    d->phi,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->D1,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->D2,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    &(d->n_lo)[1],&(d->n_hi)[1],
		    &(d->n_lo)[2],&(d->n_hi)[2],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb,&mark_D1,&mark_D2);
		    
      LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_LOCAL(d->lse_rhs,
                 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 d->phi,
		 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 d->phi0,
		 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
		 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]),
		 &use_phi0_for_sign,
		 d->index_x, d->index_y, d->index_z,
		 &(d->n_lo)[0],&(d->n_hi)[0],
		 d->narrow_band,
		 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 &mark_fb);
 
       LSM3D_TVD_RK2_STAGE1_LOCAL(d->phi_stage1,
                   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->phi,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->lse_rhs,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &dt_r,		   
		   d->index_x, d->index_y, d->index_z,
		   &(d->n_lo)[0],&(d->n_hi)[0],
		   d->narrow_band,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &mark_fb);

      /* boundary conditions */ 
      signedLinearExtrapolationBC(d->phi_stage1,g,bdry_location_idx);
      
      LSM3D_HJ_ENO2_LOCAL(d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb),&(g->khi_gb),
                    d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb),&(g->khi_gb),
		    d->phi_stage1,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->D1,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->D2,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),
		    d->index_x, d->index_y, d->index_z,
		    &(d->n_lo)[0],&(d->n_hi)[0],
		    &(d->n_lo)[1],&(d->n_hi)[1],
		    &(d->n_lo)[2],&(d->n_hi)[2],
		    d->narrow_band,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &mark_fb,&mark_D1, &mark_D2);   
 
       LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_LOCAL(d->lse_rhs,
                 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 d->phi_stage1,
		 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 d->phi0,
		 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
		 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]),
		 &use_phi0_for_sign,
		 d->index_x, d->index_y, d->index_z,
		 &(d->n_lo)[0],&(d->n_hi)[0], 
		 d->narrow_band,
                 &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		 &(g->klo_gb), &(g->khi_gb),
		 &mark_fb);
	 
       LSM3D_TVD_RK2_STAGE2_LOCAL(d->phi_next,
                   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->phi_stage1,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->phi,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->lse_rhs,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &dt_r,
		   d->index_x, d->index_y, d->index_z,
		   &(d->n_lo)[0],&(d->n_hi)[0],
		   d->narrow_band,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &mark_fb);   	 

	/* boundary conditions */ 
       signedLinearExtrapolationBC(d->phi_next,g,bdry_location_idx);
       
      /* masking enforced so that the interface stays within pore space */
       if(o->do_mask) IMPOSE_MASK_LOCAL(d->phi,d->mask,d->phi_next,g,d)	 
       else           COPY_DATA(d->phi,d->phi_next,g)
       
       t_r = t_r + dt_r;   
    }
}	 

