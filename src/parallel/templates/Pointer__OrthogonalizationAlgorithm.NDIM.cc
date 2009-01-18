/*
 * File:        Pointer__OrthogonalizationAlgorithm.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.5 $
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"

#include "OrthogonalizationAlgorithm.h"
#include "OrthogonalizationAlgorithm.cc"

template class SAMRAI::tbox::Pointer< 
  LSMLIB::OrthogonalizationAlgorithm<NDIM> 
>;
