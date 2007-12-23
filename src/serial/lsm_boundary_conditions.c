/*
 * File:        lsm_boundary_conditions.c
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date: 2006/11/01 23:22:58 $
 * Description: Implementation file for functions for imposing 
 *              boundary conditions for serial calculations
 */

#include "lsm_boundary_conditions.h"
#include "lsm_boundary_conditions2d.h"
#include "lsm_boundary_conditions3d.h"


/*============= Function definitions for boundary conditions ==============*/

void linearExtrapolationBC(
  LSMLIB_REAL *phi,
  Grid *grid,
  int bdry_location_idx)
{
  int num_dims = grid->num_dims;
  if (num_dims == 2) {

    switch (bdry_location_idx) { 
      case 0: 
      case 1: 
      case 2: 
      case 3: {
        LSM2D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &bdry_location_idx);
        break;
      }
      case 6: {
        int tmp_bdry_location_idx = 0;
        LSM2D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM2D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 7: {
        int tmp_bdry_location_idx = 2;
        LSM2D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM2D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 9: {
        int tmp_bdry_location_idx = 0;
        LSM2D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM2D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 2;
        LSM2D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM2D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      default: { /* DO NOTHING */ }

    }

  } else if (num_dims == 3) {

    switch (bdry_location_idx) { 
      case 0: 
      case 1: 
      case 2: 
      case 3: 
      case 4: 
      case 5: {
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &bdry_location_idx);
        break;
      }
      case 6: {
        int tmp_bdry_location_idx = 0;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 7: {
        int tmp_bdry_location_idx = 2;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 8: {
        int tmp_bdry_location_idx = 4;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 5;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 9: {
        int tmp_bdry_location_idx = 0;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 2;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 4;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 5;
        LSM3D_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      default: { /* DO NOTHING */ }

    }

  } /* end switch on num_dims */

}


void signedLinearExtrapolationBC(
  LSMLIB_REAL *phi,
  Grid *grid,
  int bdry_location_idx)
{
  int num_dims = grid->num_dims;
  if (num_dims == 2) {

    switch (bdry_location_idx) { 
      case 0: 
      case 1: 
      case 2: 
      case 3: {
        LSM2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &bdry_location_idx);
        break;
      }
      case 6: {
        int tmp_bdry_location_idx = 0;
        LSM2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 7: {
        int tmp_bdry_location_idx = 2;
        LSM2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 9: {
        int tmp_bdry_location_idx = 0;
        LSM2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 2;
        LSM2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      default: { /* DO NOTHING */ }

    }

  } else if (num_dims == 3) {

    switch (bdry_location_idx) { 
      case 0: 
      case 1: 
      case 2: 
      case 3: 
      case 4: 
      case 5: {
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &bdry_location_idx);
        break;
      }
      case 6: {
        int tmp_bdry_location_idx = 0;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 7: {
        int tmp_bdry_location_idx = 2;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 8: {
        int tmp_bdry_location_idx = 4;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 5;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 9: {
        int tmp_bdry_location_idx = 0;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 2;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 4;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 5;
        LSM3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      default: { /* DO NOTHING */ }

    }

  } /* end switch on num_dims */

}
 
   
void copyExtrapolationBC(
  LSMLIB_REAL *phi,
  Grid *grid,
  int bdry_location_idx)
{
  int num_dims = grid->num_dims;
  if (num_dims == 2) {

    switch (bdry_location_idx) { 
      case 0: 
      case 1: 
      case 2: 
      case 3: {
        LSM2D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &bdry_location_idx);
        break;
      }
      case 6: {
        int tmp_bdry_location_idx = 0;
        LSM2D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM2D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 7: {
        int tmp_bdry_location_idx = 2;
        LSM2D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM2D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 9: {
        int tmp_bdry_location_idx = 0;
        LSM2D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM2D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 2;
        LSM2D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM2D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      default: { /* DO NOTHING */ }

    }

  } else if (num_dims == 3) {

    switch (bdry_location_idx) { 
      case 0: 
      case 1: 
      case 2: 
      case 3: 
      case 4: 
      case 5: {
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &bdry_location_idx);
        break;
      }
      case 6: {
        int tmp_bdry_location_idx = 0;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 7: {
        int tmp_bdry_location_idx = 2;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 8: {
        int tmp_bdry_location_idx = 4;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 5;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 9: {
        int tmp_bdry_location_idx = 0;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 2;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 4;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 5;
        LSM3D_COPY_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      default: { /* DO NOTHING */ }

    }

  } /* end switch on num_dims */

}


void homogeneousNeumannBC(
  LSMLIB_REAL *phi,
  Grid *grid,
  int bdry_location_idx)
{
  copyExtrapolationBC(phi, grid, bdry_location_idx);
}
