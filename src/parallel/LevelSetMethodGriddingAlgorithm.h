/*
 * File:        LevelSetMethodGriddingAlgorithm.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.5 $
 * Modified:    $Date: 2006/01/27 16:13:20 $
 * Description: Header file for the level set method grid management class
 */
 
#ifndef included_LevelSetMethodGriddingAlgorithm_h
#define included_LevelSetMethodGriddingAlgorithm_h

/*! \class LSMLIB::LevelSetMethodGriddingAlgorithm
 *
 * \brief 
 * The LevelSetMethodGriddingAlgorithm class provides basic functionality
 * for managing grid refinement and data initialization for the the level 
 * set functions and variables involved in the computation of the velocity 
 * field.  
 *
 *
 * <h3> User-specified parameters (input database fields) </h3>
 *
 * <h4> Hierarchy Structure Input: </h4>
 *
 * - max_levels (REQUIRED)          =  integer value specifying maximum number 
 *                                     of levels allowed in the AMR patch 
 *                                     hierarchy.
 * - largest_patch_size  (REQUIRED) =  an array of integer vectors (each has 
 *                                     length = DIM) that specify the 
 *                                     dimensions of largest patch allowed 
 *                                     on each level in the hierarchy.  If
 *                                     more than max_levels entries are given,
 *                                     extra entries will be ignored. If fewer
 *                                     than max_levels entries are given, then 
 *                                     the last element in the array will be 
 *                                     used on each level without a specified 
 *                                     input value.
 * - ratio_to_coarser (REQUIRED)    =  set of (max_levels - 1) integer vectors,
 *                                     each of which indicates the ratio of 
 *                                     the index space of a patch level to 
 *                                     that of the next coarser level.  The
 *                                     input for each level must correspond to 
 *                                     the format ``level_n = vector'', where 
 *                                     n is the level number and each vector 
 *                                     must have length DIM.  
 *
 * <h4> Adaptive Refinement Input: </h4>
 * - tagging_method (OPTIONAL)      =  string array specification of the type 
 *                                     of cell-tagging used. Valid choices 
 *                                     include: ``GRADIENT_DETECTOR'' and 
 *                                     ``REFINE_BOXES''.  A combination of 
 *                                     any or all of the above may be placed
 *                                     in any order. If no input is given, no 
 *                                     tagging will be performed.  
 * - RefineBoxes (OPTIONAL)         =  input section describing the refine 
 *                                     boxes for each level.  
 *   - Level<ln>                    =  input section provides the sequence of 
 *                                     Box arrays describing where 
 *                                     user-specified refinement is to occur 
 *                                     on Level ln.
 *     - times (OPTIONAL)           =  double array specifying times at which 
 *                                     a particular box sequence is to be used.
 *     - cycles (OPTIONAL)          =  integer array specifying regrid cycles 
 *                                     at which a particular box seqence is 
 *                                     to be used.
 *     - boxes_0                    =  box array specifying refine boxes for 
 *                                     sequence 0.
 *     - boxes_1                    =  box array specifying refine boxes for 
 *                                     sequence 1.
 *     - boxes_n                    =  box array specifying refine boxes for 
 *                                     sequence n.
 *
 * <h4> Load Balancer Input: </h4>
 * - NO REQUIRED INPUT PARAMETERS (several OPTIONAL input parameters)
 *
 *
 * <h3> NOTES: </h3>
 *   - The descriptions of the input parameters was taken almost verbatim 
 *     from the class descriptions of the 
 *     SAMRAI::mesh::LevelSetMethodGriddingAlgorithm, 
 *     SAMRAI::mesh::StandardTagAndInitialize, and 
 *     SAMRAI::mesh::TagAndInitializeStrategy classes in the files
 *     LevelSetMethodGriddingAlgorithm.h, StandardTagAndInitialize.h, and
 *     TagAndInitializeStrategy.h.   For more details about the input
 *     parameters, please consult the SAMRAI documentation for these
 *     classes.
 *
 *   - For a list and description of optional gridding algorithm input 
 *     fields, see the documentation for the 
 *     SAMRAI::mesh::LevelSetMethodGriddingAlgorithm class.
 *
 *   - For a list and description of optional load balancer input 
 *     fields, see the documentation for the SAMRAI::mesh::LoadBalancer
 *     class. 
 *     
 */


#include "SAMRAI_config.h"
#include "GriddingAlgorithm.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "StandardTagAndInitialize.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "LSMLIB_config.h"
#include "LevelSetMethodGriddingStrategy.h"
#include "LevelSetFunctionIntegratorStrategy.h"
#include "LevelSetMethodVelocityFieldStrategy.h"

// SAMRAI namespaces
using namespace SAMRAI;
using namespace hier;
using namespace mesh;
using namespace tbox;


/******************************************************************
 *
 * LevelSetMethodGriddingAlgorithm Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

template<int DIM> class LevelSetMethodGriddingAlgorithm:
  public LevelSetMethodGriddingStrategy<DIM>,
  public StandardTagAndInitialize<DIM>
{
public:

  //! @{
  /*!
   ****************************************************************
   *
   * @name Constructor and destructor
   *
   ****************************************************************/

  /*!
   * This constructor for LevelSetMethodGriddingAlgorithm sets up 
   * the object to use the specified instance of the concrete subclasses 
   * of the LevelSetFunctionIntegratorStrategy and initializes the gridding 
   * algorithm to use the BergerRigoutsos box generation algorithm 
   * and the standard LoadBalancer.
   *
   * Arguments:
   *  - input_db (in):                 pointer to input database 
   *  - patch_hierarchy (in):          pointer to BasePatchHierarchy 
   *                                   object used for computation
   *  - lsm_integrator_strategy (in):  pointer to subclass of 
   *                                   LevelSetFunctionIntegratorStrategy
   *                                   that implements a gridding strategy
   *                                   based on level set function values 
   *  - object_name (in):              object name
   *
   */
  LevelSetMethodGriddingAlgorithm(
    Pointer<Database> input_db,
    Pointer< BasePatchHierarchy<DIM> > patch_hierarchy,
    Pointer< LevelSetFunctionIntegratorStrategy<DIM> > lsm_integrator_strategy,
    const string& object_name = "LevelSetMethodGriddingAlgorithm");

  /*!
   * The destructor for LevelSetMethodGriddingAlgorithm does nothing.
   */
  virtual ~LevelSetMethodGriddingAlgorithm(){};

  //! @}


  //! @{
  /*!
   ****************************************************************
   * 
   * @name Method for registering velocity field 
   *
   * Inherited from the LevelSetMethodGriddingStrategy class.
   *
   ****************************************************************/

  /*!
   * registerVelocityFieldStrategy() registers the specified instance
   * of a concrete subclass of the LevelSetMethodVelocityFieldStrategy
   * class with the LevelSetMethodGriddingAlgorithm object.  
   *
   * Arguments:     
   *  - velocity_field_strategy (in):  pointer to instance of subclass of 
   *                                   LevelSetMethodVelocityFieldStrategy 
   *                                   used to manage the variables involved 
   *                                   in the calculation of the velocity 
   *                                   field
   *
   * Return value:                     none
   *
   */
  virtual void registerVelocityFieldStrategy(
     LevelSetMethodVelocityFieldStrategy<DIM>* velocity_field_strategy);

  //! @}


  //! @{
  /*!
   ****************************************************************
   * 
   * @name Methods for managing grid configuration
   *
   * Inherited from the LevelSetMethodGriddingStrategy class.
   *
   ****************************************************************/

  /*!
   * initializePatchHierarchy() constructs the PatchHierarchy and 
   * initializes the level set functions and variables involved in 
   * the computation of the velocity field.
   *
   * Arguments:     
   *  - time (in):   simulation time that PatchHierarchy is being 
   *                 initialized
   *
   * Return value:   none
   *
   * NOTES:
   *  - all LevelSetMethodVelocityFieldStrategy objects required to 
   *    calculate the velocity field MUST be registered using 
   *    registerVelocityFieldStrategy() before invoking 
   *    initializePatchHierarchy().
   */
  virtual void initializePatchHierarchy(const double time);

  /*!
   * regridPatchHierarchy() regrids the entire PatchHierarchy and
   * reinitializes the data on the PatchHierarchy using interpolation
   * and averaging, as necessary.
   *
   * Arguments:     
   *  - time (in):   simulation time when PatchHierarchy is regrid
   *
   * Return value:   none
   *
   */
  virtual void regridPatchHierarchy(double time);

  //! @}
 
  //! @{ 
  /*!
   *******************************************************************
   * 
   * @name Methods inherited from the StandardTagAndInitialize class
   *
   *******************************************************************/

  /*!
   * initializeLevelData() initializes the data on the specifed PatchLevel
   * by invoking the initializeLevelData() method for the concrete 
   * subclasses of the LevelSetFunctionIntegratorStrategy and 
   * LevelSetMethodVelocityFieldStrategy classes that manage the time
   * evolution of the level set functions and the calculation of the 
   * velocity field, respectively.
   *
   * Arguments:     
   *  - hierarchy (in):       PatchHierarchy on which to initialize data
   *  - level_number (in):    number of PatchLevel on which initialize data
   *  - init_data_time(in):   simulation time at which data is to be
   *                          initialized
   *  - can_be_refined (in):  true if the PatchLevel can be refined; false
   *                          otherwise
   *  - initial_time (in):    true if the current time is the initial time in
   *                          the simulation; false otherwise
   *  - old_level (in):       old PatchLevel from which data is to be taken
   *                          to initialize the level_number PatchLevel
   *                          in the specified PatchHierarchy 
   *  - allocate_data (in):   true if memory for the data needs to be 
   *                          allocated; false otherwise 
   *                      
   * Return value:            none
   *                
   * NOTES:
   *  - Overrides StandardTagAndInitialize::initializeLevelData()
   *                
   */
  virtual void initializeLevelData(
    const Pointer< BasePatchHierarchy<DIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer< BasePatchLevel<DIM> > old_level = 
      Pointer< BasePatchLevel<DIM> >(NULL),
    const bool allocate_data = true);

  /*!
   * resetHierarchyConfiguration() resets the configuration of the 
   * specified PatchHierarchy by invoking the resetHierarchyConfiguration() 
   * method for the concrete subclasses of the 
   * LevelSetFunctionIntegratorStrategy and LevelSetMethodVelocityFieldStrategy
   * classes that manage the time evolution of the level set functions and 
   * the calculation of the velocity field.
   *
   * Arguments:     
   *  - hierarchy (in):       PatchHierarchy to reconfigure
   *  - coarsest_level (in):  coarsest level in hierarchy to reconfigure
   *  - finest_level (in):    finest level in hierarchy to reconfigure
   *                      
   * Return value:            none
   *                
   * NOTES:
   *  - Overrides StandardTagAndInitialize::resetHierarchyConfiguration()
   *                
   */
  void resetHierarchyConfiguration(
    const Pointer< BasePatchHierarchy<DIM> > hierarchy,
    const int coarsest_level,
    const int finest_level);

  /*!
   * tagCellsForRefinement() tags the cells on the specified PatchLevel
   * that should be refined by invoking the applyGradientDetector() 
   * and tagCellsForRefinement() methods for the concrete subclasses of 
   * the LevelSetFunctionIntegratorStrategy and the 
   * LevelSetMethodVelocityFieldStrategy classes that manage the time
   * evolution of the level set functions and the calculation of the 
   * velocity field.  In addition, the user may explicitly specify boxes 
   * which should be refined through the input file.
   *
   * Arguments:     
   *  - hierarchy (in):            PatchHierarchy which is to be regridded
   *  - level_number (in):         number of PatchLevel on which data is to be 
   *                               refined
   *  - regrid_time(in):           simulation time at which PatchLevel is to be
   *                               refined
   *  - tag_index (in):            PatchData ID containing the tag data
   *  - initial_time (in):         true if the current time is the initial time 
   *                               in the simulation; false otherwise
   *  - coarsest_sync_level (in):  true if specified PatchLevel is the coarsest
   *                               level in the hierarchy in the current
   *                               regridding process
   *  - can_be_refined (in):       true if the PatchLevel where cells are to be
   *                               tagged can be refined; false otherwise
   *  - regrid_start_time (in):    true if memory for the data needs to be 
   *                      
   * Return value:                 none
   *                
   * NOTES:
   *  - Overrides StandardTagAndInitialize::tagCellsForRefinement()
   *                
   */
  void tagCellsForRefinement(
    const Pointer< BasePatchHierarchy<DIM> > hierarchy,
    const int level_number,
    const double regrid_time,
    const int tag_index,
    const bool initial_time,
    const bool coarsest_sync_level,
    const bool can_be_refined,
    const double regrid_start_time = 0);

  //! @}


protected:

  //! @{
  /*!
   **************************************************************************
   *
   * @name Utility methods 
   *
   **************************************************************************/

  /*!
   * getFromInput() re-reads the ``tagging_method'' input from the input
   * database.  This input has already been read by the parent class,
   * but the results are not available because the relevent data members
   * are declared private in StandardTagAndInitialize.
   *
   * Arguments:     
   *  - input_db (in):  input database
   *                      
   * Return value:      none
   *                
   * NOTES:
   *  - This method does NOT override StandardTagAndInitialize::getFromInput()
   *    because getFromInput() is not declared virtual in the class
   *    StandardTagAndInitialize. 
   */
  void getFromInput(Pointer<Database> input_db);

  //! @}


  /****************************************************************
   *
   * Data Members
   *
   ****************************************************************/

  /*
   * String name for object (duplicate of 
   * StandardTagAndInitialize::d_object_name)
   */
  string d_object_name;

  /*
   * Pointer to the PatchHierarchy object
   */
  Pointer< BasePatchHierarchy<DIM> > d_patch_hierarchy;

  /*
   * Pointer to the SAMRAI::mesh::GriddingAlgorithm object
   */
  Pointer< SAMRAI::mesh::GriddingAlgorithm<DIM> > d_gridding_alg;

  /*
   * Duplicates of the booleans specifying the tagging method that
   * are data members of the StandardTagAndInitialize class.
   */
  bool d_use_refine_boxes;
  bool d_use_gradient_detector;
  bool d_use_richardson_extrapolation;

  /*
   * Pointers to the LevelSetFunctionIntegratorStrategy and 
   * LevelSetMethodVelocityFieldStrategy objects.
   */
  Pointer< LevelSetFunctionIntegratorStrategy<DIM> > d_lsm_integrator_strategy;

  Array< Pointer< LevelSetMethodVelocityFieldStrategy<DIM> > > 
    d_velocity_field_strategies;

private:
 
  /*
   * Private copy constructor to prevent use.
   * 
   * Arguments:
   *  - rhs (in):  LevelSetMethodGriddingAlgorithm object to copy
   * 
   */
  LevelSetMethodGriddingAlgorithm(
    const LevelSetMethodGriddingAlgorithm& rhs);
   
  /*
   * Private assignment operator to prevent use.
   *
   * Arguments:
   *  - rhs (in):    LevelSetMethodGriddingAlgorithm to copy
   *
   * Return value:   *this
   *
   */
  const LevelSetMethodGriddingAlgorithm& operator=(
    const LevelSetMethodGriddingAlgorithm& rhs) {
      return *this;
  }

};

} // end LSMLIB namespace

#endif
