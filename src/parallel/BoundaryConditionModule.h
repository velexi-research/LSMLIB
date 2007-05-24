/*
 * File:        BoundaryConditionModule.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.5 $
 * Modified:    $Date: 2006/10/16 20:54:25 $
 * Description: Header file for anti-periodic bc module
 */
 
#ifndef included_BoundaryConditionModule_h
#define included_BoundaryConditionModule_h

/*! \class LSMLIB::BoundaryConditionModule
 *
 * \brief 
 * The BoundaryConditionModule provides support for several common
 * boundary conditions for level set functions:  homogeneous Neumann,
 * extrapolation, and anti-periodic boundary conditions.
 *
 * In more detail, these boundary conditions are:
 * 
 * - Homogeneous Neumann boundary conditions set the ghost cells
 *   for the specified data so that the normal derivative of the 
 *   level set function at the boundary to zero.  The ghost cells
 *   are set in a manner that accounts for the type and order of 
 *   the discretization used to compute spatial derivatives.  For
 *   ENO1, ENO2, ENO3, and WENO5, this module simply fills ghost
 *   cells with data from the nearest neighbor in the interior
 *   of the computational domain. 
 *
 * - Extrapolation boundary conditions extend the level set functions
 *   into the ghost cells using simple extrapolation schemes.  There are 
 *   two types of extrapolation boundary conditions: unsigned and signed.  
 *   Unsigned extrapolation boundary conditions extend the level set 
 *   function using a simple linear extrapolation in the direction normal 
 *   to the boundary.  Signed extrapolation also extends the level set
 *   function values using a linear function, but the sign of the slope
 *   at the boundary is chosen to guarantee that the extrapolated level 
 *   set function does not have a value of zero in the extrapolation 
 *   direction.  
 * 
 * - Anti-periodic boundary conditions are periodic boundary 
 *   conditions where the sign of the function changes across the 
 *   "anti-periodic" boundary.  These boundary conditions are important 
 *   when the sign of the level set function changes across the periodic 
 *   boundary to avoid introducing artificial zero-level sets into the 
 *   level set functions at the periodic boundary.  
 *
 *   When the sign of the level set function in the interior of the 
 *   computational domain is the same on both sides of a periodic boundary, 
 *   then the values in the ghostcells outside of the physical domain are 
 *   set equal to the value of the grid cells taken from grid cells on the 
 *   other side of the computational domain.  However, when the sign of the 
 *   level set function in the interior of the computational domain changes 
 *   across a periodic boundary, the values in the ghostcells outside of the 
 *   physical domain are set to minus the value of the grid cells taken from
 *   grid cells on the other side of the computational domain.  
 *
 * <h3> NOTES: </h3>
 *
 * - In order for anti-periodic boundary conditions to be used in a 
 *   coordinate direction, the GridGeometry associated with the 
 *   PatchHierarchy MUST be set to be periodic in that coordinate 
 *   direction.  This ensures that data is properly transferred from
 *   the interior grid cells at the opposite side of the computational
 *   grid.
 *
 * - To guarantee that anti-periodic boundary conditions are correctly
 *   imposed, the level set function must truly be anti-periodic.  In
 *   particular, the initial conditions for the level set functions MUST
 *   be anti-periodic.  It is the user's responsibility to ensure that 
 *   the initial conditions for level set functions are appropriately 
 *   set when s/he chooses to impose anti-periodic boundary conditions.
 *
 */


#include "SAMRAI_config.h"
#include "BoundaryBox.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

#include "LSMLIB_config.h"
#include "LevelSetMethodToolbox.h"

// SAMRAI namespaces
using namespace SAMRAI;
using namespace hier;
using namespace tbox;


/******************************************************************
 *
 * BoundaryConditionModule Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

template<int DIM> class BoundaryConditionModule
{
public:

  /*! \enum BOUNDARY_CONDITION_TYPE
   *
   * Enumerated type for types of boundary conditions provided by
   * the BoundaryConditionModule class.
   *
   */
  typedef enum {
    NONE                        = 0,
    HOMOGENEOUS_NEUMANN         = 1,
    LINEAR_EXTRAPOLATION        = 2, 
    SIGNED_LINEAR_EXTRAPOLATION = 3,
    ANTI_PERIODIC               = 4} BOUNDARY_CONDITION_TYPE;


  //! @{
  /*!
   ****************************************************************
   *
   * @name Constructors and destructor
   *
   ****************************************************************/

  /*!
   * The standard constructor initializes the BoundaryConditionModule 
   * using the specified parameters.
   *
   * Arguments: none
   *
   * NOTES:
   *  - If the patch_hierarchy has not yet been constructed (i.e.
   *    number of levels == 0), then the BoundaryConditionModule is 
   *    set to be in an invalid state.
   */
  BoundaryConditionModule( 
    Pointer< PatchHierarchy<DIM> > patch_hierarchy,
    const IntVector<DIM>& ghostcell_width );

  /*!
   * Default constructor initializes the BoundaryConditionModule
   * to an invalid state.
   *
   * Arguments: none
   *
   */
  BoundaryConditionModule();

  /*!
   * Copy constructor. 
   *
   * Arguments:
   *  - rhs (in):  BoundaryConditionModule object to copy
   *
   */
  BoundaryConditionModule( const BoundaryConditionModule<DIM>& rhs );

  /*!
   * Destructor does nothing
   */
  virtual ~BoundaryConditionModule(){}

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for imposing boundary conditions
   *
   ****************************************************************/

  /*!
   * imposeBoundaryConditions() imposes the specified boundary conditions 
   * for phi on the entire PatchHierarchy.
   *
   * Arguments:
   *  - phi_handle (in):                PatchData handle for function on 
   *                                    which to impose boundary conditions
   *  - lower_bc (in):                  vector of integers specifying the
   *                                    type of boundary conditions to impose
   *                                    on the lower face of the computational
   *                                    domain in each coordinate direction.
   *                                    The i-th entry should contain the type
   *                                    of boundary condition to impose at 
   *                                    the lower boundary in the i-th 
   *                                    coordinate direction.
   *  - upper_bc (in):                  vector of integers specifying the
   *                                    type of boundary conditions to impose
   *                                    on the upper face of the computational
   *                                    domain in each coordinate direction.
   *                                    The i-th entry should contain the type
   *                                    of boundary condition to impose at 
   *                                    the upper boundary in the i-th 
   *                                    coordinate direction.
   *  - spatial_derivative_type (in):   type of spatial derivative 
   *                                    calculation
   *  - spatial_derivative_order (in):  order of spatial derivative
   *  - component (in):                 component of phi on which to impose 
   *                                    boundary conditions
   *                                    (default = -1)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - The spatial derivative information is only used for imposing
   *    homogeneous Neumann boundary conditions.
   *
   *  - The default behavior when no boundary conditions are supplied
   *    (i.e. either lower_bc or upper_bc is a vector of -1's) is that
   *    homogeneous Neumann boundary conditions will be imposed on all
   *    non-periodic boundaries.
   *
   *  - When component is negative, boundary conditions will be
   *    imposed on ALL of the components of phi.
   *
   */
  virtual void imposeBoundaryConditions(
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int component = -1);


  /*!
   * imposeBoundaryConditionsOnPatch() imposes the specified boundary 
   * conditions for phi on the specified patch.
   *
   * Arguments:
   *  - patch (in):                     Patch on which set boundary conditions
   *  - phi_handle (in):                PatchData handle for function on 
   *                                    which to impose boundary conditions
   *  - lower_bc (in):                  vector of integers specifying the
   *                                    type of boundary conditions to impose
   *                                    on the lower face of the computational
   *                                    domain in each coordinate direction.
   *                                    The i-th entry should contain the type
   *                                    of boundary condition to impose at 
   *                                    the lower boundary in the i-th 
   *                                    coordinate direction.
   *  - upper_bc (in):                  vector of integers specifying the
   *                                    type of boundary conditions to impose
   *                                    on the upper face of the computational
   *                                    domain in each coordinate direction.
   *                                    The i-th entry should contain the type
   *                                    of boundary condition to impose at 
   *                                    the upper boundary in the i-th 
   *                                    coordinate direction.
   *  - spatial_derivative_type (in):   type of spatial derivative 
   *                                    calculation
   *  - spatial_derivative_order (in):  order of spatial derivative
   *  - component (in):                 component of phi on which to impose 
   *                                    boundary conditions
   *                                    (default = -1)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - The spatial derivative information is only used for imposing
   *    homogeneous Neumann boundary conditions.
   *
   *  - The default behavior when no boundary conditions are supplied
   *    (i.e. either lower_bc or upper_bc is a vector of -1's) is that
   *    homogeneous Neumann boundary conditions will be imposed on all
   *    non-periodic boundaries.
   *
   *  - When component is negative, boundary conditions will be
   *    imposed on ALL of the components of phi.
   *
   */
  virtual void imposeBoundaryConditionsOnPatch(
    Patch<DIM>& patch,
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int component = -1);


  /*!
   * imposeAntiPeriodicBCs() imposes "anti-periodic" boundary 
   * conditions at periodic boundaries where the sign of a
   * level set function changes.
   *
   * Arguments:
   *  - phi_handle (in):  PatchData handle for function on which 
   *                      to impose anti-periodic boundary 
   *                      conditions
   *  - lower_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the lower face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the lower boundary in the i-th 
   *                      coordinate direction.
   *  - upper_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the upper face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the upper boundary in the i-th 
   *                      coordinate direction.
   *  - component (in):   component of phi on which to impose 
   *                      anti-periodic boundary conditions
   *                      (default = -1)
   *
   * Return value:        none
   *
   * NOTES:
   *  - Anti-periodic boundary conditions are only imposed for those
   *    directions that are specified by lower_bc and upper_bc AND that
   *    are periodic directions for the GridGeometry object associated 
   *    with the PatchHierarchy set in the constructor or by 
   *    resetHierarchyConfiguration().  If a direction is specified
   *    to be anti-periodic by the lower_bc and upper_bc variables but 
   *    is not a periodic direction for the GridGeometry object, then 
   *    that direction is NOT treated as an anti-periodic direction.
   *
   *  - This method assumes that ghostcells for ALL patches have already
   *    been correctly filled for true periodic boundary conditions (e.g.
   *    via SAMRAI communication routines such as fillData()).  This method 
   *    merely corrects the sign of the ghostcell values when the level set
   *    function is anti-periodic across the periodic boundary.
   *  
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - When component is negative, anti-periodic boundary conditions will 
   *    be imposed on ALL of the components of phi.
   *
   */
  virtual void imposeAntiPeriodicBCs(
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const int component = -1);


  /*!
   * imposeAntiPeriodicBCsOnPatch() imposes "anti-periodic" boundary 
   * conditions on the specified Patch at periodic boundaries where 
   * the sign of a level set function changes.
   *
   * Arguments:
   *  - patch (in):       Patch on which set boundary conditions
   *  - phi_handle (in):  PatchData handle for function on which 
   *                      to impose anti-periodic boundary 
   *                      conditions
   *  - lower_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the lower face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the lower boundary in the i-th 
   *                      coordinate direction.
   *  - upper_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the upper face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the upper boundary in the i-th 
   *                      coordinate direction.
   *  - component (in):   component of phi on which to impose 
   *                      anti-periodic boundary conditions
   *                      (default = -1)
   *
   * Return value:        none
   *
   * NOTES:
   *  - Anti-periodic boundary conditions are only imposed for those
   *    directions that are specified by lower_bc and upper_bc AND that
   *    are periodic directions for the GridGeometry object associated 
   *    with the PatchHierarchy set in the constructor or by 
   *    resetHierarchyConfiguration().  If a direction is specified
   *    to be anti-periodic by the lower_bc and upper_bc variables but is
   *    not a periodic direction for the GridGeometry object, then 
   *    that direction is NOT treated as an anti-periodic direction.
   *
   *  - This method assumes that ghostcells for ALL patches have already
   *    been correctly filled for true periodic boundary conditions (e.g.
   *    via SAMRAI communication routines such as fillData()).  This method 
   *    merely corrects the sign of the ghostcell values when the level set
   *    function is anti-periodic across the periodic boundary.
   *  
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - When component is negative, anti-periodic boundary conditions will 
   *    be imposed on ALL of the components of phi.
   *
   */
  virtual void imposeAntiPeriodicBCsOnPatch(
    Patch<DIM>& patch,
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const int component = -1);


  /*!
   * imposeHomogeneousNeumannBCs() imposes homogeneous Neumann boundary 
   * conditions at the specified boundary locations.
   *
   * Arguments:
   *  - phi_handle (in):                PatchData handle for function on 
   *                                    which to impose homogeneous Neumann 
   *                                    boundary conditions
   *  - lower_bc (in):                  vector of integers specifying the
   *                                    type of boundary conditions to impose
   *                                    on the lower face of the computational
   *                                    domain in each coordinate direction.
   *                                    The i-th entry should contain the type
   *                                    of boundary condition to impose at 
   *                                    the lower boundary in the i-th 
   *                                    coordinate direction.
   *  - upper_bc (in):                  vector of integers specifying the
   *                                    type of boundary conditions to impose
   *                                    on the upper face of the computational
   *                                    domain in each coordinate direction.
   *                                    The i-th entry should contain the type
   *                                    of boundary condition to impose at 
   *                                    the upper boundary in the i-th 
   *                                    coordinate direction.
   *  - spatial_derivative_type (in):   spatial derivative type
   *  - spatial_derivative_order (in):  spatial derivative order
   *  - component (in):                 component of phi on which to 
   *                                    impose homogeneous Neumann 
   *                                    boundary conditions
   *                                    (default = -1)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - It is the user's responsibility to ensure that the data 
   *    associated with phi_handle is sufficient for imposing a
   *    homogeneous Neumann boundary condition for the specified
   *    discretization of the spatial derivative.
   *
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - When component is negative, boundary conditions will be imposed 
   *    for ALL components of phi.
   *
   */
  virtual void imposeHomogeneousNeumannBCs(
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int component = -1);


  /*!
   * imposeHomogeneousNeumannBCsOnPatch() imposes homogeneous Neumann 
   * boundary conditions at the specified boundary locations on the
   * specified Patch.
   *
   * Arguments:
   *  - patch (in):                     Patch on which set boundary conditions
   *  - phi_handle (in):                PatchData handle for function on 
   *                                    which to impose homogeneous Neumann 
   *                                    boundary conditions
   *  - lower_bc (in):                  vector of integers specifying the
   *                                    type of boundary conditions to impose
   *                                    on the lower face of the computational
   *                                    domain in each coordinate direction.
   *                                    The i-th entry should contain the type
   *                                    of boundary condition to impose at 
   *                                    the lower boundary in the i-th 
   *                                    coordinate direction.
   *  - upper_bc (in):                  vector of integers specifying the
   *                                    type of boundary conditions to impose
   *                                    on the upper face of the computational
   *                                    domain in each coordinate direction.
   *                                    The i-th entry should contain the type
   *                                    of boundary condition to impose at 
   *                                    the upper boundary in the i-th 
   *                                    coordinate direction.
   *  - spatial_derivative_type (in):   spatial derivative type
   *  - spatial_derivative_order (in):  spatial derivative order
   *  - component (in):                 component of phi on which to 
   *                                    impose homogeneous Neumann 
   *                                    boundary conditions
   *                                    (default = -1)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - It is the user's responsibility to ensure that the data 
   *    associated with phi_handle is sufficient for imposing a
   *    homogeneous Neumann boundary condition for the specified
   *    discretization of the spatial derivative.
   *
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - When component is negative, boundary conditions will be imposed 
   *    for ALL components of phi.
   *
   */
  virtual void imposeHomogeneousNeumannBCsOnPatch(
    Patch<DIM>& patch,
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int component = -1);


  /*!
   * imposeLinearExtrapolationBCs() imposes linear extrapolation
   * boundary conditions at the specified boundary locations.
   *
   * Arguments:
   *  - phi_handle (in):  PatchData handle for function on which to
   *                      impose linear extrapolation boundary 
   *                      conditions
   *  - lower_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the lower face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the lower boundary in the i-th 
   *                      coordinate direction.
   *  - upper_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the upper face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the upper boundary in the i-th 
   *                      coordinate direction.
   *  - component (in):   component of phi on which to impose linear
   *                      extrapolation boundary conditions
   *                      (default = -1)
   *
   * Return value:        none
   *
   * NOTES:
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - When component is negative, boundary conditions will be imposed 
   *    for ALL components of phi.
   *
   */
  virtual void imposeLinearExtrapolationBCs(
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const int component = -1);

  /*!
   * imposeLinearExtrapolationBCsOnPatch() imposes linear extrapolation
   * boundary conditions at the specified boundary locations on the
   * specified Patch.
   *
   * Arguments:
   *  - patch (in):       Patch on which set boundary conditions
   *  - phi_handle (in):  PatchData handle for function on which to
   *                      impose linear extrapolation boundary 
   *                      conditions
   *  - lower_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the lower face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the lower boundary in the i-th 
   *                      coordinate direction.
   *  - upper_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the upper face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the upper boundary in the i-th 
   *                      coordinate direction.
   *  - component (in):   component of phi on which to impose linear
   *                      extrapolation boundary conditions
   *                      (default = -1)
   *
   * Return value:        none
   *
   * NOTES:
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - When component is negative, boundary conditions will be imposed 
   *    for ALL components of phi.
   *
   */
  virtual void imposeLinearExtrapolationBCsOnPatch(
    Patch<DIM>& patch,
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const int component = -1);

  /*!
   * imposeSignedLinearExtrapolationBCs() imposes signed-linear 
   * extrapolation boundary conditions at the specified boundary 
   * locations.
   *
   * Arguments:
   *  - phi_handle (in):  PatchData handle for function on which to
   *                      impose signed linear extrapolation boundary 
   *                      conditions
   *  - lower_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the lower face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the lower boundary in the i-th 
   *                      coordinate direction.
   *  - upper_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the upper face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the upper boundary in the i-th 
   *                      coordinate direction.
   *  - component (in):   component of phi on which to impose signed 
   *                      linear extrapolation boundary conditions
   *                      (default = -1)
   *
   * Return value:        none
   *
   * NOTES:
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - When component is negative, boundary conditions will be imposed 
   *    for ALL components of phi.
   *
   */
  virtual void imposeSignedLinearExtrapolationBCs(
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const int component = -1);

  /*!
   * imposeSignedLinearExtrapolationBCsOnPatch() imposes signed-linear 
   * extrapolation boundary conditions at the specified boundary 
   * locations on the specified Patch.
   *
   * Arguments:
   *  - patch (in):       Patch on which set boundary conditions
   *  - phi_handle (in):  PatchData handle for function on which to
   *                      impose signed linear extrapolation boundary 
   *                      conditions
   *  - lower_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the lower face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the lower boundary in the i-th 
   *                      coordinate direction.
   *  - upper_bc (in):    vector of integers specifying the
   *                      type of boundary conditions to impose
   *                      on the upper face of the computational
   *                      domain in each coordinate direction.
   *                      The i-th entry should contain the type
   *                      of boundary condition to impose at 
   *                      the upper boundary in the i-th 
   *                      coordinate direction.
   *  - component (in):   component of phi on which to impose signed 
   *                      linear extrapolation boundary conditions
   *                      (default = -1)
   *
   * Return value:        none
   *
   * NOTES:
   *  - The number of ghostcells for the data associated with phi_handle
   *    MUST be equal to the ghostcell_width in the specified in the 
   *    constructor or resetHierarchyConfiguration() method. It is the 
   *    user's responsibility to ensure that the current configuration
   *    of the BoundaryConditionModule is compatible with the phi data.
   *
   *  - When component is negative, boundary conditions will be imposed 
   *    for ALL components of phi.
   *
   */
  virtual void imposeSignedLinearExtrapolationBCsOnPatch(
    Patch<DIM>& patch,
    const int phi_handle,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const int component = -1);

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for managing the grid configuration
   *
   ****************************************************************/

  /*!
   * resetHierarchyConfiguration() resets the configuration of the
   * BoundaryConditionModule to be consistent with the specified 
   * PatchHierarchy.  In particular,  it computes ALL of the boundary 
   * boxes (both periodic and non-periodic) that need to be filled.  
   * The SAMRAI library internally carries out the same calculation, 
   * but it is necessary to repeat this calculation in order to 
   * impose anti-periodic boundary conditions at periodic boundaries 
   * across which level set functions change sign.
   *
   * Arguments:
   *  - patch_hierarchy (in):        PatchHierarchy to reconfigure
   *  - coarsest_level (in):         coarsest level in hierarchy to reconfigure
   *  - finest_level (in):           finest level in hierarchy to reconfigure
   *  - ghostcell_width (in):        width of ghostcell required for
   *                                 boundary conditions
   *
   * Return value:                   none
   *
   * NOTES:
   *  - If the patch_hierarchy has not yet been constructed (i.e.
   *    number of levels == 0), then the BoundaryConditionModule is 
   *    set to be in an invalid state.
   *
   */
  virtual void resetHierarchyConfiguration(
    const Pointer< PatchHierarchy<DIM> > patch_hierarchy,
    const int coarsest_level,
    const int finest_level,
    const IntVector<DIM>& ghostcell_width);

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Operators
   *
   ****************************************************************/

  /*!
   * Assignment operator. 
   *
   * Arguments:
   *  - rhs (in):    BoundaryConditionModule to copy
   *
   * Return value:   *this
   *
   */
  virtual inline const BoundaryConditionModule& operator=(
    const BoundaryConditionModule<DIM>& rhs)
  {
    d_patch_hierarchy = rhs.d_patch_hierarchy;
    d_ghostcell_width = rhs.d_ghostcell_width;
    d_geom_periodic_dirs = rhs.d_geom_periodic_dirs;
    d_boundary_boxes = rhs.d_boundary_boxes;
    d_touches_boundary = rhs.d_touches_boundary;
    return *this;
  }

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Static helper methods
   *
   ****************************************************************/

  /*!
   * computeIndexSpaceOfNearestGhostLayer() computes the index space
   * of the layer of ghost cells nearest the interior of the computational
   * domain.
   * 
   * Arguments:
   *  - nearest_ghost_layer_lo (out):  lower corner of box that represents
   *                                   the index space of ghost-cell layer
   *                                   nearest the boundary specified
   *                                   box to fill 
   *  - nearest_ghost_layer_hi (out):  upper corner of box that represents
   *                                   the index space of ghost-cell layer
   *                                   nearest the boundary specified
   *                                   box to fill 
   *  - bdry_type (in):                boundary type (see SAMRAI 
   *                                   documentation for definition)
   *  - bdry_location_idx (in):        boundary location index (see SAMRAI
   *                                   documentation for definition)
   *  - fillbox (in):                  box (of ghost-cells) to fill
   *  
   * Return value:                     none
   * 
   */
  static void computeIndexSpaceOfNearestGhostLayer(
    IntVector<DIM>& nearest_ghost_layer_lo,
    IntVector<DIM>& nearest_ghost_layer_hi,
    const int bdry_type,
    const int bdry_location_idx,
    const Box<DIM>& fillbox);
    
  /*!
   * computeIndexOffset() computes the offset in the data array between
   * the ghost cell nearest the interior of the computational domain
   * and the first interior cell.
   *
   * Arguments:
   *  - bdry_type (in):          boundary type (see SAMRAI documentation
   *                             for definition)
   *  - bdry_location_idx (in):  boundary location index (see SAMRAI
   *                             documentation for definition)
   *  - ghostbox (in):           box of ghost-cells to fill
   *
   * Return value:               none
   *
   */
  static int computeIndexOffset(
    const int bdry_type,
    const int bdry_location_idx,
    const Box<DIM>& ghostbox);
  
  //! @}


protected:

  /****************************************************************
   *
   * Data members
   *
   ****************************************************************/

  // pointer to PatchHierarchy
  Pointer< PatchHierarchy<DIM> > d_patch_hierarchy;

  // parameters for imposing anti-periodic BCs
  IntVector<DIM> d_ghostcell_width;
  IntVector<DIM> d_geom_periodic_dirs;
  Array< Array< Array< BoundaryBox<DIM> > > > d_boundary_boxes;
  Array< Array<bool> > d_touches_boundary;
  
};

} // end LSMLIB namespace 

#endif
