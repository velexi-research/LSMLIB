/*
 * File:        LevelSetMethodVelocityFieldStrategy.NDIM.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"

#include "LevelSetMethodVelocityFieldStrategy.h"
#include "LevelSetMethodVelocityFieldStrategy.cc"

template class LSMLIB::LevelSetMethodVelocityFieldStrategy<NDIM>;
