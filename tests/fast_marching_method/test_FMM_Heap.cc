/*
 * Test program for FMM_Heap
 *
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE. This file is part of the LSMLIB package. It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution. No part of the LSMLIB
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

#include <stdlib.h>                 // for rand, NULL, srand, RAND_MAX
#include <time.h>                   // for time

#include "gtest/gtest-message.h"    // for Message
#include "gtest/gtest-test-part.h"  // for TestPartResult
#include "gtest/gtest_pred_impl.h"  // for Test, ASSERT_EQ, ASSERT_LE, ...

#include "lsmlib_config.h"
#include "FMM_Heap.h"

/*
 * Constants
 */
#define TEST_DIM (2)

/*
 * Tests
 */
TEST(FMMHeapTest, FMMHeap)
{
    int i,j;
    int grid_idx[FMM_HEAP_MAX_NDIM];

    /*
     * Heap Test code
     */
    FMM_Heap *fmm_heap = FMM_Heap_createHeap(0,0,0);
    ASSERT_TRUE(FMM_Heap_isEmpty(fmm_heap));
    EXPECT_EQ(FMM_Heap_getHeapSize(fmm_heap), 0);

    int N = 4;  // number of cells in each grid direction
    srand((int)time(0));  // seed random number generator

    // insert some nodes
    for (i = 0; i<N; i++) {
        for (j = 0; j<N; j++) {
            LSMLIB_REAL value;
            int node_handle;

            grid_idx[0] = i;
            grid_idx[1] = j;
            value = 1.0*rand()/RAND_MAX;
            node_handle = FMM_Heap_insertNode(fmm_heap,grid_idx,value);
        }
    }

    // extract some nodes
    LSMLIB_REAL prev_val = -1;
    for (i = 0; i < 2*N; i++) {
        FMM_HeapNode moved_node;
        int moved_handle;
        FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,&moved_node,&moved_handle);

        // check results
        ASSERT_LE(prev_val, root.value);
        ASSERT_EQ(root.heap_pos, 0);

        // update prev_val
        prev_val = root.value;
    }

    // insert some more nodes
    for (i = 0; i<N; i++) {
        for (j = 0; j<N; j++) {
            LSMLIB_REAL value = 1.0*rand()/RAND_MAX;
            int node_handle = FMM_Heap_insertNode(fmm_heap,grid_idx,value);

            grid_idx[0] = i+N;
            grid_idx[1] = j+N;
        }
    }

    // extract some more nodes without the moved_node argument
    prev_val = -1;  // reset prev_val
    for (i = 0; i < N; i++) {
        int moved_handle;
        FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,NULL,&moved_handle);

        // check results
        ASSERT_LE(prev_val, root.value);
        ASSERT_EQ(root.heap_pos, 0);

        // update prev_val
        prev_val = root.value;
    }

    // extract some more nodes without the moved_handle argument
    prev_val = -1;  // reset prev_val
    for (i = 0; i < N; i++) {
        FMM_HeapNode moved_node;
        FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,&moved_node,NULL);

        // check results
        ASSERT_LE(prev_val, root.value);
        ASSERT_EQ(root.heap_pos, 0);

        // update prev_val
        prev_val = root.value;
    }

    // extract some more nodes without the moved_node and moved_handle arguments
    prev_val = -1;  // reset prev_val
    for (i = 0; i < N; i++) {
        FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,NULL,NULL);

        // check results
        ASSERT_LE(prev_val, root.value);
        ASSERT_EQ(root.heap_pos, 0);

        // update prev_val
        prev_val = root.value;
    }

    // update some nodes
    for (i = 0; i<N; i++) {
        LSMLIB_REAL value = i+1;
        FMM_HeapNode node = FMM_Heap_getNode(fmm_heap,i);

        FMM_Heap_updateNode(fmm_heap,i,value);
        node = FMM_Heap_getNode(fmm_heap,i);

        // check results
        ASSERT_EQ(node.value, i+1);
    }
    for (i = N; i<2*N; i++) {
        LSMLIB_REAL value = -i;
        FMM_HeapNode node = FMM_Heap_getNode(fmm_heap,i);

        FMM_Heap_updateNode(fmm_heap,i,value);
        node = FMM_Heap_getNode(fmm_heap,i);

        // check results
        ASSERT_EQ(node.value, -i);
    }

    // extract remaining nodes
    prev_val = -1e100;  // reset prev_val
    while (FMM_Heap_getHeapSize(fmm_heap) > 0) {
        FMM_HeapNode moved_node;
        int moved_handle;
        FMM_HeapNode root = FMM_Heap_extractMin(fmm_heap,&moved_node,&moved_handle);

        // check results
        ASSERT_LE(prev_val, root.value);
        ASSERT_EQ(root.heap_pos, 0);

        // update prev_val
        prev_val = root.value;
    }

    // clean up memory
    FMM_Heap_destroyHeap(fmm_heap);
}
