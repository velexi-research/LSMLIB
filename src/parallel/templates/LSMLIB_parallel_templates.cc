/*
 * File:        Array__Array__Array__bool.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Array.C"
#include "FieldExtensionAlgorithm.h"
#include "FieldExtensionAlgorithm.cc"
#include "boost/shared_ptr.hpp"
#include "LevelSetMethodAlgorithm.h"
#include "LevelSetMethodAlgorithm.cc"
#include "boost/config.hpp"
#include "LevelSetMethodVelocityFieldStrategy.h"
#include "LevelSetMethodVelocityFieldStrategy.cc"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/RefineSchedule.C"

template class SAMRAI::tbox::Array<boost::shared_ptr<LSMLIB::FieldExtensionAlgorithm> >;

template class SAMRAI::tbox::Array<boost::shared_ptr<LSMLIB::LevelSetMethodVelocityFieldStrategy> >;

template class SAMRAI::tbox::Array<SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> > >;

template class boost::shared_ptr<LSMLIB::LevelSetMethodAlgorithm>; 
