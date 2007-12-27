/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* LSMLIB headers */
#include "LSMLIB_config.h"
#include "lsm_level_set_evolution3d.h"
#include "lsm_spatial_derivatives3d.h"
#include "lsm_utilities3d.h"
#include "lsm_tvd_runge_kutta3d.h"
#include "lsm_reinitialization3d.h"
#include "lsm_initialization3d.h"
#include "lsm_geometry3d.h"
#include "lsm_boundary_conditions3d.h"

/* LSMLIB Serial package headers */
#include "lsm_boundary_conditions.h"
#include "lsm_macros.h"

/* Local headers */
#include "curvature_model_top.h"
#include "curvature_model3d.h"

#define DT_MIN_TO_CORRECT 1e-5
#define DT_MIN            0.001


void curvatureModelMedium3dMainLoop(
     Options          *options,
     LSM_DataArrays   *data_arrays,
     Grid             *grid,
     FILE             *fp_out)
{
  LSMLIB_REAL   cfl_number = 0.5;
  
  /* time variables */
  LSMLIB_REAL   t, dt, dt_sub, max_H, dt_corr;
  LSMLIB_REAL   tplot, dt_min, dt_max;
  LSMLIB_REAL   tmax_r = 5*grid->dx[0]; /* max time for reinitialization */
  
  LSMLIB_REAL   max_abs_err, eps, eps_stop;
 
  LSMLIB_REAL   zero = 0.0;
  LSMLIB_REAL   vel_n, vol_phi, vol_max, vol_phi_prev, rel_vol_diff;
  int      nx, nxy;  
  
  int      bdry_location_idx = 9; /* extrapolate all boundaries */
  
  int      OUTER_STEP, INNER_STEP, TOTAL_STEP;
  int      reinit_steps, last_reinit_step, ave_reinit_steps;
  
  /* writing shortcuts */
  Grid             *g = grid;
  LSM_DataArrays   *d = data_arrays;
  Options          *o = options;
  
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
   		
  OUTER_STEP = 0; INNER_STEP = 0; TOTAL_STEP = 0;
  last_reinit_step = 0;
  reinit_steps = ave_reinit_steps = 0;
  			 
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
      
      SET_DATA_TO_CONSTANT(d->lse_rhs,g,zero)    
     
      if(o->a > 0)
      {  
         /* Compute upwinding gradient approximations */
          LSM3D_HJ_ENO2(d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
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
		    &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		    &(g->klo_fb), &(g->khi_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));   
         vel_n = o->a;
	 
	 LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS(d->lse_rhs,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &vel_n,
		    &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		    &(g->klo_fb), &(g->khi_fb));
                  
	 /* figure out dt for hyperbolic term */
	 LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT(&dt,&vel_n,
		    d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_minus, d->phi_y_minus,  d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		    &(g->klo_fb), &(g->khi_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]),		    
		    &cfl_number);  
      }
      else dt = tplot;
     
      if( o->b > 0)
      {
	/* Compute derivatives needed for curvature term*/
	LSM3D_CENTRAL_GRAD_ORDER2(d->phi_x, d->phi_y, d->phi_z,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(g->ilo_D1_fb), &(g->ihi_D1_fb), 
		    &(g->jlo_D1_fb), &(g->jhi_D1_fb),
		    &(g->klo_D1_fb), &(g->khi_D1_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));	
	LSM3D_CENTRAL_GRAD_ORDER2(d->phi_xx, d->phi_xy, d->phi_xz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(g->ilo_D2_fb), &(g->ihi_D2_fb), 
		    &(g->jlo_D2_fb), &(g->jhi_D2_fb),
		    &(g->klo_D2_fb), &(g->khi_D2_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));	
	LSM3D_CENTRAL_GRAD_ORDER2(d->phi_xy, d->phi_yy, d->phi_yz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_y,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(g->ilo_D2_fb), &(g->ihi_D2_fb), 
		    &(g->jlo_D2_fb), &(g->jhi_D2_fb),
		    &(g->klo_D2_fb), &(g->khi_D2_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
	LSM3D_CENTRAL_GRAD_ORDER2(d->phi_xz, d->phi_yz, d->phi_zz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_z,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(g->ilo_D2_fb), &(g->ihi_D2_fb), 
		    &(g->jlo_D2_fb), &(g->jhi_D2_fb),
		    &(g->klo_D2_fb), &(g->khi_D2_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));		    
	
	LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS(d->lse_rhs,
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
		    &(g->ilo_D2_fb), &(g->ihi_D2_fb), 
		    &(g->jlo_D2_fb), &(g->jhi_D2_fb),
		    &(g->klo_D2_fb), &(g->khi_D2_fb));
	 
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
      
      LSM3D_TVD_RK2_STAGE1(d->phi_stage1,
                   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->phi,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->lse_rhs,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		   &(g->klo_fb), &(g->khi_fb),
		   &dt);
      /* boundary conditions */	   
      signedLinearExtrapolationBC(d->phi_stage1,g,bdry_location_idx);
      
      /* masking enforced so that the interface stays within pore space */
      if(o->do_mask) IMPOSE_MASK(d->phi_stage1,d->mask,d->phi_stage1,g);   	  

      SET_DATA_TO_CONSTANT(d->lse_rhs,g,zero)
      
      if(o->a)
      {
	  LSM3D_HJ_ENO2(d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
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
		    &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		    &(g->klo_fb), &(g->khi_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
	 LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS(d->lse_rhs,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x_minus, d->phi_y_minus, d->phi_z_minus,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &vel_n,
		    &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		    &(g->klo_fb), &(g->khi_fb));	    
      }
      
      if( o->b )
      {
	LSM3D_CENTRAL_GRAD_ORDER2(d->phi_x,d->phi_y,d->phi_z,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_stage1,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(g->ilo_D1_fb), &(g->ihi_D1_fb), 
		    &(g->jlo_D1_fb), &(g->jhi_D1_fb),
		    &(g->klo_D1_fb), &(g->khi_D1_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));	
	LSM3D_CENTRAL_GRAD_ORDER2(d->phi_xx,d->phi_xy,d->phi_xz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_x,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(g->ilo_D2_fb), &(g->ihi_D2_fb), 
		    &(g->jlo_D2_fb), &(g->jhi_D2_fb),
		    &(g->klo_D2_fb), &(g->khi_D2_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));	
	LSM3D_CENTRAL_GRAD_ORDER2(d->phi_xy,d->phi_yy,d->phi_yz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_y,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(g->ilo_D2_fb), &(g->ihi_D2_fb), 
		    &(g->jlo_D2_fb), &(g->jhi_D2_fb),
		    &(g->klo_D2_fb), &(g->khi_D2_fb), 
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
	LSM3D_CENTRAL_GRAD_ORDER2(d->phi_xz, d->phi_yz, d->phi_zz,
                    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    d->phi_z,
		    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		    &(g->klo_gb), &(g->khi_gb),
		    &(g->ilo_D2_fb), &(g->ihi_D2_fb), 
		    &(g->jlo_D2_fb), &(g->jhi_D2_fb),
		    &(g->klo_D2_fb), &(g->khi_D2_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));		    
	//Compute_and_add_curvature3d_lse_rhs(g,p,curv_tmp,grad_mag2)
	LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS(d->lse_rhs,
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
		    &(g->ilo_D2_fb), &(g->ihi_D2_fb), 
		    &(g->jlo_D2_fb), &(g->jhi_D2_fb),
		    &(g->klo_D2_fb), &(g->khi_D2_fb));
      }
     
      LSM3D_TVD_RK2_STAGE2(d->phi_next,
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
		   &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		   &(g->klo_fb), &(g->khi_fb),
		   &dt);

      /* boundary conditions */
      signedLinearExtrapolationBC(d->phi_next,g,bdry_location_idx);	 
	 
      /* masking enforced so that the interface stays within pore space */	 
       if(o->do_mask) IMPOSE_MASK(d->phi,d->mask,d->phi_next,g)  
       else           COPY_DATA(d->phi,d->phi_next,g)
       
      dt_sub = dt_sub + dt;
   } /*inner loop */
  
   t = t + dt_sub;
    
   if(o->do_reinit)
   {  /* periodic reinitialization of the level set function */      
      reinitializeMedium3d(d,g,o,tmax_r);
      
      ave_reinit_steps += TOTAL_STEP - last_reinit_step;
      last_reinit_step = TOTAL_STEP;
      reinit_steps++;   	 
   }
    
   /* compute max abs error */
   LSM3D_MAX_NORM_DIFF(&max_abs_err,d->phi,
            &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	    &(g->klo_gb), &(g->khi_gb),
	    d->phi_prev,
	    &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	    &(g->klo_gb), &(g->khi_gb),
	    &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
	    &(g->klo_fb), &(g->khi_fb));
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
  } /* outer loop */
  

  /* Print out statistics on reinitialization */ 
  ave_reinit_steps =ceil( (LSMLIB_REAL)(ave_reinit_steps) / (LSMLIB_REAL)(reinit_steps) );
  fprintf(fp_out,"\nTotal steps %d  Reinit. steps %d",TOTAL_STEP,reinit_steps);
  fprintf(fp_out,"\nReinitialized on average every %d steps.\n",ave_reinit_steps);
  
}


/* 
*  reinitializeMedium3d() reinitializes the level set function using the second
*  order accuracy ENO and TVD RK routines.
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
   

void reinitializeMedium3d(
     LSM_DataArrays *data_arrays,
     Grid           *grid,
     Options        *options,
     LSMLIB_REAL         tmax_r)
{  
    LSMLIB_REAL cfl_number = 0.5;
    LSMLIB_REAL t_r, dt_r;
    
    int    use_phi0_for_sign = 0;    
    int    bdry_location_idx = 9; /* all boundaries */
   
      /* writing shortcuts */
    Grid             *g = grid;
    LSM_DataArrays   *d = data_arrays;
    Options          *o = options;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
      
    COPY_DATA(d->phi0,d->phi,g)
    
    while(t_r < tmax_r )
    {
      LSM3D_HJ_ENO2(d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
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
		    &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		    &(g->klo_fb), &(g->khi_fb),
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));   
		    		    
      LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS(d->lse_rhs,
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
		 &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		 &(g->klo_fb), &(g->khi_fb),
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]),
		 &use_phi0_for_sign);
 
       LSM3D_TVD_RK2_STAGE1(d->phi_stage1,
                   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->phi,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   d->lse_rhs,
		   &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
		   &(g->klo_gb), &(g->khi_gb),
		   &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		   &(g->klo_fb), &(g->khi_fb),
		   &dt_r);
      /* boundary conditions */
      signedLinearExtrapolationBC(d->phi_stage1,g,bdry_location_idx);	 	 	    
      
      LSM3D_HJ_ENO2(d->phi_x_plus, d->phi_y_plus, d->phi_z_plus,
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
		    &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		    &(g->klo_fb), &(g->khi_fb), 
		    &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
       LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS(d->lse_rhs,
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
		 &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		 &(g->klo_fb), &(g->khi_fb),
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]),
		 &use_phi0_for_sign);
	 
       LSM3D_TVD_RK2_STAGE2(d->phi_next,
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
		   &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
		   &(g->klo_fb), &(g->khi_fb),
		   &dt_r);
   	   
       /* boundary conditions */
       signedLinearExtrapolationBC(d->phi_next,g,bdry_location_idx); 
      
       /* masking enforced so that the interface stays within pore space */
       if(o->do_mask) IMPOSE_MASK(d->phi,d->mask,d->phi_next,g)	 
       else           COPY_DATA(d->phi,d->phi_next,g)
       
       t_r = t_r + dt_r;   
    }
}	 

