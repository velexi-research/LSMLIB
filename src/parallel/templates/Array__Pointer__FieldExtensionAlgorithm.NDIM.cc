/*
 * File:        Array__Pointer__FieldExtensionAlgorithm.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "FieldExtensionAlgorithm.h"
#include "tbox/Array.h"
#include "tbox/Array.C"
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"

namespace SAMRAI {

template class tbox::Array< tbox::Pointer< 
  LSMLIB::FieldExtensionAlgorithm<NDIM> 
> >;

}
