/*
 * File:        Array__Array__Array__BoundaryBox.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "BoundaryBox.h"
#include "BoundaryBox.C"
#include "tbox/Array.h"
#include "tbox/Array.C"

namespace SAMRAI {

template class tbox::Array< tbox::Array< 
  tbox::Array< hier::BoundaryBox<NDIM> > 
> >;

}
