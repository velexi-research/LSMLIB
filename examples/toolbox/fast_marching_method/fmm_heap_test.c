/*
 * File:        fmm_heap_test.c
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.2 $
 * Modified:    $Date: 2006/12/04 13:54:23 $
 * Description: test program for the FMM_Heap "object"
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "FMM_Heap.h"

/************************************************************************
 * Helper function declarations
 ************************************************************************/
#define TEST_DIM (2)
void printGridIndex(int grid_idx[FMM_HEAP_MAX_NDIM]);


/************************************************************************
 *                                                                    
 * Driver code for testing the FMM_Heap class.                       
 *                                                                  
 ************************************************************************/

int main( int argc, char *argv[])
{
  int i,j;
  int grid_idx[FMM_HEAP_MAX_NDIM];

  /* 
   * Heap Test code
   */
  FMM_Heap *fmm_heap = FMM_Heap_createHeap(0,0,0);
  FMM_Heap_printHeapData(fmm_heap);

  int N = 4;  // number of cells in each grid direction
  srand((int)time(0));  // seed random number generator

  // insert some nodes
  printf("Inserting some nodes...\n");
  for (i = 0; i<N; i++) {
    for (j = 0; j<N; j++) {
      double value;
      int node_handle;

      grid_idx[0] = i;
      grid_idx[1] = j;
      value = 1.0*rand()/RAND_MAX;
      node_handle = FMM_Heap_insertNode(fmm_heap,grid_idx,value);
      printf("Node handle = %d, ", node_handle);
      printf("Grid Index = "); printGridIndex(grid_idx); printf(", ");
      printf("Value = %g, \n", value); 
      printf("   Heap Size = %d, ", FMM_Heap_getHeapSize(fmm_heap));
      printf("Heap Mem Size = %d\n", FMM_Heap_getHeapMemSize(fmm_heap));
    } 
  } 

  // extract some nodes
  double prev_val = -1;
  printf("\nExtracting some nodes...\n");
  for (i = 0; i < 2*N; i++) {
    FMM_HeapNode moved_node; 
    int moved_handle;
    FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,&moved_node,&moved_handle);
    printf("---------------------\n");
    printf("Root: \n");
    printf("Grid Index = "); printGridIndex(root.grid_idx); printf(", ");
    printf("Value = %g, ", root.value); 
    printf("Heap Position = %d\n", root.heap_pos);
    printf("Moved Node: \n");
    printf("Grid Index = "); printGridIndex(moved_node.grid_idx); printf(", ");
    printf("Value = %g, ", moved_node.value);
    printf("Heap Position = %d\n", moved_node.heap_pos); 
    printf("   Node Handle = %d\n", moved_handle); 
    printf("Heap Size = %d, ", FMM_Heap_getHeapSize(fmm_heap));
    printf("Heap Mem Size = %d\n", FMM_Heap_getHeapMemSize(fmm_heap));
    if (prev_val > root.value) 
      printf("ERROR!!!  Heap Property Failed!!!\n");
    if (0 != root.heap_pos) 
      printf("ERROR!!!  heap_pos field set incorrectly!!!\n");
    printf("---------------------\n");

    // update prev_val
    prev_val = root.value;

  } 

  // insert some more nodes
  printf("\nInserting some more nodes...\n");
  for (i = 0; i<N; i++) {
    for (j = 0; j<N; j++) {
      double value = 1.0*rand()/RAND_MAX;
      int node_handle = FMM_Heap_insertNode(fmm_heap,grid_idx,value);

      grid_idx[0] = i+N;
      grid_idx[1] = j+N;

      printf("Node handle = %d, ", node_handle);
      printf("Grid Index = "); printGridIndex(grid_idx); printf(", ");
      printf("Value = %g\n", value); 
      printf("   Heap Size = %d, ", FMM_Heap_getHeapSize(fmm_heap));
      printf("Heap Mem Size = %d\n", FMM_Heap_getHeapMemSize(fmm_heap));
    } 
  }

  // extract some more nodes without the moved_node argument
  printf("\nExtracting some more nodes...no moved_node argument\n");
  prev_val = -1;  // reset prev_val
  for (i = 0; i < N; i++) {
    int moved_handle;
    FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,NULL,&moved_handle);
    printf("---------------------\n");
    printf("Root: \n");
    printf("Grid Index = "); printGridIndex(root.grid_idx); printf(", ");
    printf("Value = %g, ", root.value); 
    printf("Heap Position = %d\n", root.heap_pos); 
    printf("Moved Node: \n");
    printf("   Node Handle = %d\n", moved_handle); 
    printf("Heap Size = %d, ", FMM_Heap_getHeapSize(fmm_heap));
    printf("Heap Mem Size = %d\n", FMM_Heap_getHeapMemSize(fmm_heap));
    if (prev_val > root.value) 
      printf("ERROR!!!  Heap Property Failed!!!\n");
    if (0 != root.heap_pos) 
      printf("ERROR!!!  heap_pos field set incorrectly!!!\n");
    printf("---------------------\n");

    // update prev_val
    prev_val = root.value;
  }

  // extract some more nodes without the moved_handle argument
  printf("\nExtracting some more nodes...no moved_handle argument\n"); 
  prev_val = -1;  // reset prev_val
  for (i = 0; i < N; i++) {
    FMM_HeapNode moved_node; 
    FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,&moved_node,NULL);
    printf("---------------------\n");
    printf("Root: \n");
    printf("Grid Index = "); printGridIndex(root.grid_idx); printf(", ");
    printf("Value = %g, ", root.value); 
    printf("Heap Position = %d\n", root.heap_pos); 
    printf("Moved Node: \n");
    printf("Grid Index = "); printGridIndex(moved_node.grid_idx); printf(", ");
    printf("Value = %g, ", moved_node.value); 
    printf("Heap Position = %d\n", moved_node.heap_pos); 
    printf("Heap Size = %d, ", FMM_Heap_getHeapSize(fmm_heap));
    printf("Heap Mem Size = %d\n", FMM_Heap_getHeapMemSize(fmm_heap));
    if (prev_val > root.value) 
      printf("ERROR!!!  Heap Property Failed!!!\n");
    if (0 != root.heap_pos) 
      printf("ERROR!!!  heap_pos field set incorrectly!!!\n");
    printf("---------------------\n");

    // update prev_val
    prev_val = root.value;
  }

  // extract some more nodes without the moved_node and moved_handle arguments
  printf("\nExtracting some more nodes...no moved_node\n");
  printf(" AND moved_handle arguments\n");
  prev_val = -1;  // reset prev_val
  for (i = 0; i < N; i++) {
    FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,NULL,NULL);
    printf("---------------------\n");
    printf("Root: \n");
    printf("Grid Index = "); printGridIndex(root.grid_idx); printf(", ");
    printf("Value = %g, ", root.value);
    printf("Heap Position = %d\n", root.heap_pos);
    printf("Heap Size = %d, ", FMM_Heap_getHeapSize(fmm_heap));
    printf("Heap Mem Size = %d\n ", FMM_Heap_getHeapMemSize(fmm_heap));
    if (prev_val > root.value) 
      printf("ERROR!!!  Heap Property Failed!!!\n");
    if (0 != root.heap_pos) 
      printf("ERROR!!!  heap_pos field set incorrectly!!!\n");
    printf("---------------------\n");

    // update prev_val
    prev_val = root.value;
  }
 
  // update some nodes
  printf("\nUpdating value of some nodes...\n");
  for (i = 0; i<N; i++) {
    double value = i+1;
    FMM_HeapNode node = FMM_Heap_getNode(fmm_heap,i);
    printf("Before:  \n");
    printf("   Node handle = %d, ", i);
    printf("Grid Index = "); printGridIndex(node.grid_idx); printf(", ");
    printf("Value = %g, ", node.value);
    printf("Heap Position = %d\n", node.heap_pos); 
    FMM_Heap_updateNode(fmm_heap,i,value);
    node = FMM_Heap_getNode(fmm_heap,i);
    printf("After:  \n");
    printf("   Node handle = %d, ", i); 
    printf("Grid Index = "); printGridIndex(node.grid_idx); printf(", ");
    printf("Value = %g, ", node.value); 
    printf("Heap Position = %d\n", node.heap_pos); 
    printf("   Heap Size = %d, ", FMM_Heap_getHeapSize(fmm_heap));
    printf("Heap Mem Size = %d\n", FMM_Heap_getHeapMemSize(fmm_heap));

    if (node.value != i+1) 
      printf("ERROR: FMM_Heap_updateNode() failed!!!\n");
  }
  for (i = N; i<2*N; i++) {
    double value = -i;
    FMM_HeapNode node = FMM_Heap_getNode(fmm_heap,i);
    printf("Before:  \n");
    printf("   Node handle = %d, ", i);
    printf("Grid Index = "); printGridIndex(node.grid_idx); printf(", ");
    printf("Value = %g, ", node.value);
    printf("Heap Position = %d\n", node.heap_pos);
    FMM_Heap_updateNode(fmm_heap,i,value);
    node = FMM_Heap_getNode(fmm_heap,i);
    printf("After:  \n");
    printf("   Node handle = %d, ", i); 
    printf("Grid Index = "); printGridIndex(node.grid_idx); printf(", ");
    printf("Value = %g, ", node.value); 
    printf("Heap Position = %d\n", node.heap_pos); 
    printf("   Heap Size = %d, ", FMM_Heap_getHeapSize(fmm_heap));
    printf("Heap Mem Size = %d\n", FMM_Heap_getHeapMemSize(fmm_heap));

    if (node.value != -i) 
      printf("ERROR: FMM_Heap_updateNode() failed!!!\n");
  }

  // extract remaining nodes
  printf("\nExtracting remaining nodes...\n");
  prev_val = -1e100;  // reset prev_val
  while (FMM_Heap_getHeapSize(fmm_heap) > 0) {
    FMM_HeapNode moved_node; 
    int moved_handle;
    FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,&moved_node,&moved_handle);
    printf("---------------------\n");
    printf("Root: \n");
    printf("Grid Index = "); printGridIndex(root.grid_idx); printf(", ");
    printf("Value = %g, ", root.value); 
    printf("Heap Position = %d\n", root.heap_pos); 
    printf("Moved Node: \n");
    printf("Grid Index = "); printGridIndex(moved_node.grid_idx); printf(", ");
    printf("Value = %g, ",  moved_node.value); 
    printf("Heap Position = %d\n", moved_node.heap_pos);
    printf("   Node Handle = %d\n", moved_handle);
    printf("Heap Size = %d, ", FMM_Heap_getHeapSize(fmm_heap));
    printf("Heap Mem Size = %d\n", FMM_Heap_getHeapMemSize(fmm_heap));
    if (prev_val > root.value) 
      printf("ERROR!!!  Heap Property Failed!!!\n");
    if (0 != root.heap_pos) 
      printf("ERROR!!!  heap_pos field set incorrectly!!!\n");
    printf("---------------------\n");

    // update prev_val
    prev_val = root.value;
  }
 
  // clean up memory
  FMM_Heap_destroyHeap(fmm_heap); 

  return(0);
}

void printGridIndex(int grid_idx[FMM_HEAP_MAX_NDIM])
{
  int i;

  printf("(");
  for (i = 0; i<TEST_DIM; i++) {
    printf("%d", grid_idx[i]);
    if (i<TEST_DIM-1) printf(",");
  } 
  printf(")");
}
