/*
 * File:        Array__Array__Array__bool.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "tbox/Array.h"
#include "tbox/Array.C"

template class SAMRAI::tbox::Array< 
  SAMRAI::tbox::Array< SAMRAI::tbox::Array<bool> >
>;
