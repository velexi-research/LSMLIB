/*
 * File:        Array__Array__Pointer__RefineSchedule.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "RefineSchedule.h"
#include "tbox/Array.h"
#include "tbox/Array.C"
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"

namespace SAMRAI {

template class tbox::Array< tbox::Array< 
  tbox::Pointer< xfer::RefineSchedule<NDIM> > 
> >;

}
