/*
 * File:        demo_grid_management.c
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.3 $
 * Modified:    $Date: 2006/12/05 21:16:17 $
 * Description: Test code for grid management functions.
 */

#include <stdio.h>
#include "LSMLIB_config.h"
#include "lsm_file.h"
#include "lsm_grid.h"

int main(void)
{
  Grid *g_original, *g_from_ascii_file, *g_from_binary_file;
  int num_dims;
  int grid_dims[3];
  LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy_type;
  LSMLIB_REAL x_lo[3];
  LSMLIB_REAL x_hi[3];


  printf("\n************ Testing 2D Grid management functions ************\n\n");

  num_dims = 2;
  grid_dims[0] = 10; grid_dims[1] = 20;
  accuracy_type = VERY_HIGH;
  x_lo[0] = -1; x_lo[1] = -1;
  x_hi[0] = 1; x_hi[1] = 1;

  g_original = createGridSetGridDims(num_dims, grid_dims, 
                                     x_lo, x_hi, accuracy_type);
  printf("====================== Original Grid =========================\n");
  printGrid(g_original,stdout); 
  printf("==============================================================\n");

  writeGridToAsciiFile(g_original,"demo_grid_2d.ascii",NO_ZIP);
  g_from_ascii_file = readGridFromAsciiFile("demo_grid_2d.ascii");
  printf("=================== Grid From ASCII File =====================\n");
  printGrid(g_from_ascii_file,stdout);
  printf("==============================================================\n");

  writeGridToBinaryFile(g_original,"demo_grid_2d.binary",GZIP);
  g_from_binary_file = readGridFromBinaryFile("demo_grid_2d.binary");
  printf("=================== Grid From Binary File ====================\n");
  printGrid(g_from_binary_file,stdout);
  printf("==============================================================\n");

  printf("\n**************** End 2D Grid management demo ****************\n\n");


  printf("\n************ Testing 3D Grid management functions ************\n\n");

  num_dims = 3;
  grid_dims[0] = 10; grid_dims[1] = 20; grid_dims[2] = 30;
  accuracy_type = VERY_HIGH;
  x_lo[0] = -1; x_lo[1] = -1; x_lo[2] = -1;
  x_hi[0] = 1; x_hi[1] = 1; x_hi[2] = 1;

  g_original = createGridSetGridDims(num_dims, grid_dims, 
                                     x_lo, x_hi, accuracy_type);
  printf("====================== Original Grid =========================\n");
  printGrid(g_original,stdout); 
  printf("==============================================================\n");

  writeGridToAsciiFile(g_original,"demo_grid_3d.ascii",NO_ZIP);
  g_from_ascii_file = readGridFromAsciiFile("demo_grid_3d.ascii");
  printf("=================== Grid From ASCII File =====================\n");
  printGrid(g_from_ascii_file,stdout);
  printf("==============================================================\n");

  writeGridToBinaryFile(g_original,"demo_grid_3d.binary",GZIP);
  g_from_binary_file = readGridFromBinaryFile("demo_grid_3d.binary");
  printf("=================== Grid From Binary File ====================\n");
  printGrid(g_from_binary_file,stdout);
  printf("==============================================================\n");

  printf("\n**************** End 3D Grid management demo ****************\n\n");

  return 0;
}
