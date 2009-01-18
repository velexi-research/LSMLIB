/*
 * File:        ReinitializationAlgorithm.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"

#include "ReinitializationAlgorithm.h"
#include "ReinitializationAlgorithm.cc"

template class LSMLIB::ReinitializationAlgorithm<NDIM>;
