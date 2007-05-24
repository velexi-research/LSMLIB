/*
 * File:        Pointer__BoundaryConditionModule.NDIM.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/10/04 19:11:39 $
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"

#include "BoundaryConditionModule.h"
#include "BoundaryConditionModule.cc"

template class SAMRAI::tbox::Pointer< LSMLIB::BoundaryConditionModule<NDIM> >;
