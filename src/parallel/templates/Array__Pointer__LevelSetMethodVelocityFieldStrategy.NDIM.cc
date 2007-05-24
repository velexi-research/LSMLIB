/*
 * File:        Array__Pointer__LevelSetMethodVelocityFieldStrategy.NDIM.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.5 $
 * Modified:    $Date: 2006/10/02 00:47:34 $
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "LevelSetMethodVelocityFieldStrategy.h"
#include "LevelSetMethodVelocityFieldStrategy.cc"
#include "tbox/Array.h"
#include "tbox/Array.C"
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"

namespace SAMRAI {

template class tbox::Array< tbox::Pointer< 
  LSMLIB::LevelSetMethodVelocityFieldStrategy<NDIM> 
> >;

}
