/*
 * File:        Array__Array__Array__bool.cc
 * Copyright:   (c) 2005 Kevin T. Chu
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/10/04 19:11:39 $
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "tbox/Array.h"
#include "tbox/Array.C"

template class SAMRAI::tbox::Array< 
  SAMRAI::tbox::Array< SAMRAI::tbox::Array<bool> >
>;
