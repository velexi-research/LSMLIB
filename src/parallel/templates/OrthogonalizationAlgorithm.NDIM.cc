/*
 * File:        OrthogonalizationAlgorithm.NDIM.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date: 2006/10/02 00:47:35 $
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"

#include "OrthogonalizationAlgorithm.h"
#include "OrthogonalizationAlgorithm.cc"

template class LSMLIB::OrthogonalizationAlgorithm<NDIM>;
