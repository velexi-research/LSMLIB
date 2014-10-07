/*
 * File:        Pointer__LevelSetMethodGriddingAlgorithm.NDIM.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
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
