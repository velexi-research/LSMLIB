/*
 * File:        lsm_reinitialization_medium2d.c
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Support file for reinitialization.
 */
 
/* System headers */
#include <stdio.h>
#include <stdlib.h>

/* Headers for LSM evolution eqns */
#include "lsm_level_set_evolution2d.h"
#include "lsm_spatial_derivatives2d.h"
#include "lsm_math_utils2d.h"
#include "lsm_tvd_runge_kutta2d.h"
#include "lsm_reinitialization2d.h"
#include "lsm_geometry2d.h"

#include "lsm_macros.h"
#include "lsm_boundary_conditions.h"
#include "lsm_reinitialization_medium2d.h"

/* Main loop for reinitialization. */
     
void lsm2dReinitializationMedium(
  LSM_DataArrays *lsm_arrays,
  Grid *grid,
  Options *options)
{   
  LSMLIB_REAL cfl_number = 0.5;
  LSMLIB_REAL tmax_r, t_r, dt_r;

  int    use_phi0_for_sign = 1;
  int    idx;

  /* shorten writing */
  LSM_DataArrays *l = lsm_arrays;
  Grid           *g = grid;
  Options        *o = options;
  
  tmax_r = o->tmax;

  t_r = 0;
  dt_r = cfl_number * (g->dx)[0];

  COPY_DATA(l->phi0,l->phi,g);

  while(t_r < tmax_r )
  { 
    LSM2D_HJ_ENO2(l->phi_x_plus, l->phi_y_plus,
             &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
             l->phi_x_minus, l->phi_y_minus,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->D1,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->D2,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
	     &((g->dx)[0]),
	     &((g->dx)[1]));

    LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS(l->lse_rhs,
             &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi0,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi_x_plus, l->phi_y_plus,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi_x_minus, l->phi_y_minus,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
	     &((g->dx)[0]), &((g->dx)[1]),
	     &use_phi0_for_sign);

    LSM2D_TVD_RK2_STAGE1(l->phi_stage1,
             &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->lse_rhs,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
	     &dt_r);

    /* boundary conditions */	  
    signedLinearExtrapolationBC(l->phi_stage1, g, ALL_BOUNDARIES);

    LSM2D_HJ_ENO2(l->phi_x_plus, l->phi_y_plus,
             &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
             l->phi_x_minus, l->phi_y_minus,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi_stage1,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->D1,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->D2,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
	     &((g->dx)[0]),
	     &((g->dx)[1]));	 
    LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS(l->lse_rhs,
             &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi_stage1,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi0,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi_x_plus, l->phi_y_plus,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi_x_minus, l->phi_y_minus,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
	     &((g->dx)[0]), &((g->dx)[1]),
	     &use_phi0_for_sign);

    LSM2D_TVD_RK2_STAGE2(l->phi_next,
             &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi_stage1,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->phi,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     l->lse_rhs,
	     &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	     &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
	     &dt_r);

    /* boundary conditions */	   
    signedLinearExtrapolationBC(l->phi_next, g, ALL_BOUNDARIES);

     /* update l->phi */ 
    if( o->do_mask ) 
    {
      IMPOSE_MASK(l->phi,l->mask,l->phi_next,g)
    }  
    else
    {
      COPY_DATA(l->phi,l->phi_next,g)
    }
    t_r = t_r + dt_r;   
  }
}	 

