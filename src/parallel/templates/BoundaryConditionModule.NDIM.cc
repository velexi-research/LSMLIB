/*
 * File:        BoundaryConditionModule.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/10/04 19:11:39 $
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"

#include "BoundaryConditionModule.h"
#include "BoundaryConditionModule.cc"

template class LSMLIB::BoundaryConditionModule<NDIM>;
