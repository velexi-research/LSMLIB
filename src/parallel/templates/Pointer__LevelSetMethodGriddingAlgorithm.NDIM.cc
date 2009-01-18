/*
 * File:        Pointer__LevelSetMethodGriddingAlgorithm.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"

#include "LevelSetMethodGriddingAlgorithm.h"
#include "LevelSetMethodGriddingAlgorithm.cc"

template class SAMRAI::tbox::Pointer< 
  LSMLIB::LevelSetMethodGriddingAlgorithm<NDIM> 
>;
