/*
 * File:        Pointer__FieldExtensionAlgorithm.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.5 $
 * Modified:    $Date: 2006/10/02 00:47:35 $
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"

#include "FieldExtensionAlgorithm.h"
#include "FieldExtensionAlgorithm.cc"

template class SAMRAI::tbox::Pointer< LSMLIB::FieldExtensionAlgorithm<NDIM> >;
