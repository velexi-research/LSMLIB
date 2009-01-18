/*
 * File:        OrthogonalizationAlgorithm.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"

#include "OrthogonalizationAlgorithm.h"
#include "OrthogonalizationAlgorithm.cc"

template class LSMLIB::OrthogonalizationAlgorithm<NDIM>;
