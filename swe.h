//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/examples/swe/swe.h $
// Package:     SAMRAI applications
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Numerical routines for single patch in swe equation ex.
//
 
#ifndef included_sweXD
#define included_sweXD

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "SAMRAI/algs/HyperbolicPatchStrategy.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Array.h"

//#include "tbox/AbstractStream.h"

#include <string>
using namespace std;
#define included_String

/**
 * The swe class provides routines for a sample application code that
 * solves the shallow water equations with rotation. The code illustrates the 
 * manner in which a code employing the standard Berger/Oliger AMR algorithm 
 * for explicit hydrodynamics can be used in the SAMRAI framework.
 * This class is derived from the algs::HyperbolicPatchStrategy abstract base 
 * class which defines the bulk of the interface between the hyperbolic 
 * intergration algorithm provided by SAMRAI and the numerical routines
 * specific to swe.  In particular, this class provides routines which
 * maybe applied to any patch in an AMR patch hierarchy. 
 * 
 * The numerical routines model the shallow water equations using  
 * explicit timestepping and a second order unsplit Godunov method.  
 * The primary numerical quantities forming the state vector are 
 * the conserved quantities depth and depth*velocity (mass flux)
 */

using namespace SAMRAI;

class swe : 
   public tbox::Serializable,
   public algs::HyperbolicPatchStrategy,
   public appu::BoundaryUtilityStrategy,
   public appu::VisDerivedDataStrategy
{
public:
   //==================================================================
   // constructor:
   //    - set default params
   //    - register swe object for restart
   //    - restart or get variables from input
   //==================================================================
   
   swe(
		const string& object_name,
		const tbox::Dimension& dim,
		tbox::Pointer<tbox::Database> input_db,
		tbox::Pointer<geom::CartesianGridGeometry > grid_geom);  

   //==================================================================
   // destructor
   //==================================================================
     
   ~swe();
 
   ///
   ///  The following routines:
   ///
   ///      registerModelVariables(),
   ///      initializeDataOnPatch(),
   ///      computeStableDtOnPatch(),
   ///      computeFluxesOnPatch(),
   ///      conservativeDifferenceOnPatch(),
   ///      tagGradientDetectorCells(),
   ///      tagRichardsonExtrapolationCells()
   ///
   ///  are concrete implementations of functions declared in the
   ///  algs::HyperbolicPatchStrategy abstract base class.
   ///

   /**
    * Register swe model variables with algs::HyperbolicLevelIntegrator
    * according to variable registration function provided by the integrator.
    * In other words, variables are registered according to their role
    * in the integration process (e.g., time-dependent, flux, etc.).
    * This routine also registers variables for plotting with the
    * Vis writer (Vizamrai or VisIt).  
    */
   //==================================================================
   //        registerModelVariables
   //
   //    - register model fields with the Hyperbolic Level Integrator
   //      variables are registered according to their role  
   //      e.g.:   TIME_DEP or FLUX
   //      also specify intrinsic or user defined routines for 
   //      coarsening and refinement of fields
   //   
   //    - register variables for output
   //==================================================================
   void registerModelVariables(algs::HyperbolicLevelIntegrator* integrator);

   //==================================================================
   //        setupLoadBalancer     
   //
   // called by gridding, this routine allows users to specify
   // a specific non-uniform load balance.  Here, from the Euler example
   // code, a weight of 1 is applied everywhere so that the result
   // should be identical to a uniform load
   //
   //==================================================================
   void setupLoadBalancer(algs::HyperbolicLevelIntegrator* integrator,
                          mesh::GriddingAlgorithm* gridding_algorithm); 

//======================================================================
//    set initial state on a patch (no restart) 
//       note: during sim, new patches get data through interpolation
//             this routine is initial model state only
//
//    this is a concrete implementation of a function from
//    the algs::HyperbolicPatchStrategy abstract base class
//
//    the actual implementation in swe.C calls a fortran function 
//    to perform the actual initialization
//======================================================================

   void initializeDataOnPatch(hier::Patch& patch,
                              const double data_time,
                              const bool initial_time);

//======================================================================
//         computeStableDtOnPatch
//
//     - compute stable time step on a patch 
//     - returns stable DT (double)
//     - I believe there is an internally defined CFL that is modifiable
//       through a separate call to the hyperbolic integrator
//
//     this is a concrete implementation of a function from
//     the algs::HyperbolicPatchStrategy abstract base class
//
//     the actual implementation in swe.C calls a fortran function 
//     to perform the actual initialization
//======================================================================

   double computeStableDtOnPatch(hier::Patch& patch,
                                 const bool initial_time,
                                 const double dt_time);

//======================================================================
//          computeFluxesOnPatch
//
//    - compute the fluxes on each cell face on a patch
//
//     this is a concrete implementation of a function from
//     the algs::HyperbolicPatchStrategy abstract base class
//
//     the actual implementation in swe.C calls a fortran function 
//     to perform the actual initialization
//======================================================================

   void computeFluxesOnPatch(hier::Patch& patch, 
                             const double time, 
                             const double dt);

//======================================================================
//         conservativeDifferenceOnPatch
//
//    update solution variables to new time level using flux difference
//    fields are maintained in primitive form, but update is made
//    in consdervative form.  This is the way it was done in the Euler
//    example, but is to me a bit strange.  Either way you are reconstructing 
//    the different forms at some point.
//    
//    this is a concrete implementation of a function from  
//    the algs::HyperbolicPatchStrategy abstract base class

//     the actual implementation in swe.C calls a fortran function 
//     to perform the flux difference.  This function first converts 
//     to conservative form, summs the fluxes (actually the difference)
//     then converts back to primitive
//======================================================================

   void conservativeDifferenceOnPatch(hier::Patch& patch,
                                      const double time,
                                      const double dt,
                                      bool at_syncronization);

   /**
    * Tag cells for refinement using gradient detector.
    */
   void tagGradientDetectorCells(hier::Patch& patch,
      const double regrid_time,
      const bool initial_error,
      const int tag_indx,
      const bool uses_richardson_extrapolation_too);

   /**
    * Tag cells for refinement using Richardson extrapolation.
    */
   void tagRichardsonExtrapolationCells(
      hier::Patch& patch,
      const int error_level_number,
      const tbox::Pointer<hier::VariableContext> coarsened_fine,
      const tbox::Pointer<hier::VariableContext> advanced_coarse,
      const double regrid_time,
      const double deltat,
      const int error_coarsen_ratio,
      const bool initial_error,
      const int tag_index,
      const bool uses_gradient_detector_too);

   ///
   ///  The following routines:   
   ///
   ///      setPhysicalBoundaryConditions(),
   ///      getRefineOpStencilWidth(),
   ///      postprocessRefine()
   ///
   ///  are concrete implementations of functions declared in the
   ///  RefinePatchStrategy abstract base class.
   ///

   /**
    * Set the data in ghost cells corresponding to physical boundary
    * conditions.  Specific boundary conditions are determined by 
    * information specified in input file and numerical routines.
    */
   void setPhysicalBoundaryConditions(hier::Patch& patch,
                                      const double fill_time,
                                      const hier::IntVector& 
                                      ghost_width_to_fill);

   /**
    * Return stencil width of conservative linear interpolation operations.
    */
   //hier::IntVector getRefineOpStencilWidth() const;


   /**
    * Write state of swe object to the given database for restart.
    * 
    * This routine is a concrete implementation of the function 
    * declared in the tbox::Serializable abstract base class.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
    * boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face (in 3D) or edge (in 2D) to which the boundary 
    * condition applies.
    */
   void readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                       string& db_name,
                                       int bdry_location_index);

   /**
    * Register a VisIt data writer so this class will write
    * plot files that may be postprocessed with the VisIt 
    * visualization tool.
    */
#ifdef HAVE_HDF5
   void registerVisItDataWriter( 
      tbox::Pointer<appu::VisItDataWriter > viz_writer);
#endif

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class appu::VisDerivedDataStrategy.  It computes derived
    * plot quantities registered with the Vizamrai or VisIt data 
    * writers from data  that is maintained on each patch in the 
    * hierarchy.  In particular, it writes the plot quantity 
    * identified by the string variable name to the specified 
    * double buffer on the patch in the given region.  The depth_id
    * integer argument indicates which entry in the "depth" of the    
    * vector is being written; for a scalar quantity, this may be
    * ignored.  For a vector quantity, it may be used to compute
    * the quantity at the particular depth (e.g. mom[depth_id] = 
    * rho * vel[depth_id]).  The boolean return value specifies
    * whether or not derived data exists on the patch.  Generally,
    * this will be TRUE.  If the derived data does NOT exist on 
    * the patch, return FALSE.  
    */
   bool packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const string& variable_name,
      int zeta_id) const;

   ///
   ///  The following routines are specific to the swe class and
   ///  are not declared in any base class.
   ///

   /**
    * Reset physical boundary values in special cases, such as when
    * using symmetric (i.e., reflective) boundary conditions.
    */
   void boundaryReset(hier::Patch& patch,
                      pdat::FaceData<double>& traced_left,
                      pdat::FaceData<double>& traced_right) const;

   /**
    * Print all data members for swe class.
    */
   void printClassData(ostream& os) const;

   /*
    * Dump data in intersection of 1-dimensional "pencil box" to file
    * with given name.  The direction corresponds to the axis of the
    * pencil box in the domain.  Data dumped by this routine is 
    * readable by Matlab.
    */
   void writeData1dPencil(const tbox::Pointer<hier::Patch > patch,
                          const hier::Box& pencil_box,
                          const int idir,
                          ostream& file);

private:
   /*
    * These private member functions read data from input and restart.
    * When beginning a run from a restart file, all data members are read
    * from the restart file.   If the boolean flag "is_from_restart"
    * is true when reading from input, some restart values may be
    * overridden by those in the input file.
    *
    * An assertion results if the database pointer is null.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db,
                     bool is_from_restart);
   void getFromRestart();

   void readStateDataEntry(tbox::Pointer<tbox::Database> db,
                           const string& db_name,
                           int array_indx,
                           tbox::Array<double>& depth,
                           tbox::Array<double>& bathy,
                           tbox::Array<double>& veldepth);

   /*
    * Private member function to check correctness of boundary data.
    */
   void checkBoundaryData(int btype,
                          const hier::Patch& patch,
                          const hier::IntVector& ghost_width_to_fill,
                          const tbox::Array<int>& scalar_bconds,
                          const tbox::Array<int>& vector_bconds) const;

   /*
    * The object name is used for error/warning reporting and also as a 
    * string label for restart database entries. 
    */
   string d_object_name;

   /*
    * We cache pointers to the grid geometry and Vis data writers
    * to set up initial data, set physical boundary conditions,
    * and register plot variables.  We also cache a pointer to the 
    * plot context passed to the variable registration routine.
    */
   tbox::Pointer<geom::CartesianGridGeometry > d_grid_geometry;
#ifdef HAVE_HDF5
   tbox::Pointer<appu::VisItDataWriter > d_visit_writer;
#endif
   tbox::Pointer<hier::VariableContext> d_plot_context; 

	/*
 	* Problem spatial dimension.
 	*/
	const tbox::Dimension d_dim;

   /*
    * Data items used for nonuniform load balance, if used.
    */
   tbox::Pointer< pdat::CellVariable<double> > d_workload_variable;
   int d_workload_data_id;
   bool d_use_nonuniform_workload;

   //swe state vectors containing time-dep sol [d,u,v]
   tbox::Pointer< pdat::CellVariable<double> > d_depth;
   tbox::Pointer< pdat::CellVariable<double> > d_veldepth;
   tbox::Pointer< pdat::CellVariable<double> > d_scalar;

   //auxiliary vectors => bathymetry, drycell, addmass
   tbox::Pointer< pdat::CellVariable<double> > d_bathy;
   tbox::Pointer< pdat::CellVariable<double> > d_drycell;
   tbox::Pointer< pdat::CellVariable<double> > d_addvol;

   //Pointer to the left and right running fluctuations
   tbox::Pointer< pdat::FaceVariable<double> > d_fluxm;
   tbox::Pointer< pdat::FaceVariable<double> > d_fluxp;

   //physical constants
   double d_rho;
   double d_gravity;

 
   //store number of ghosts cells for state vars and for fluxes
   hier::IntVector d_nghosts;
   hier::IntVector d_fluxghosts;

   /*
    * Indicator for problem type and initial conditions
    */
   string d_data_problem;
   int d_data_problem_int;

   /*
    * Input for SPHERE problem
    */ 
   /*
   double d_radius;
   double d_center[NDIM];
   double d_depth_inside; 
   double d_veldepth_inside[NDIM]; 
   double d_zepth_outside; 
   double d_veldepth_outside[NDIM]; 
   */

   /*
    * Input for PIECEWISE_CONSTANT_*  and STEP problems
    */
   /*
   int d_number_of_intervals;
   tbox::Array<double> d_front_position;
   tbox::Array<double> d_interval_depth; 
   tbox::Array<double> d_interval_veldepth; 
   */

   /*
    * Boundary condition cases and boundary values.
    * Options are: FLOW, REFLECT, DIRICHLET
    * and variants for nodes and edges.
    *
    * Input file values are read into these arrays.
    */
   tbox::Array<int> d_master_bdry_edge_conds;
   tbox::Array<int> d_master_bdry_node_conds;

   /*
    * Boundary condition cases for scalar and vector (i.e., depth > 1)
    * variables.  These are post-processed input values and are passed
    * to the boundary routines.
    */
   tbox::Array<int> d_scalar_bdry_edge_conds;
   tbox::Array<int> d_vector_bdry_edge_conds;

   tbox::Array<int> d_scalar_bdry_node_conds;
   tbox::Array<int> d_vector_bdry_node_conds;

   tbox::Array<int> d_node_bdry_edge;


   /*
    * Arrays of face (3d) or edge (2d) boundary values for DIRICHLET case.
    */
   tbox::Array<double> d_bdry_edge_depth;
   tbox::Array<double> d_bdry_edge_bathy;
   tbox::Array<double> d_bdry_edge_veldepth;

   /*
    * Refinement criteria parameters for gradient detector and
    * Richardson extrapolation. 
    */
   tbox::Array<string> d_refinement_criteria;
   tbox::Array<double> d_depth_grad_tol;
   tbox::Array<double> d_depth_grad_time_max;
   tbox::Array<double> d_depth_grad_time_min;
   tbox::Array<double> d_depth_rich_tol;
   tbox::Array<double> d_depth_rich_time_max;
   tbox::Array<double> d_depth_rich_time_min;

   tbox::Array<double> d_bathy_grad_tol;
   tbox::Array<double> d_bathy_grad_time_max;
   tbox::Array<double> d_bathy_grad_time_min;
   tbox::Array<double> d_bathy_rich_tol;
   tbox::Array<double> d_bathy_rich_time_max;
   tbox::Array<double> d_bathy_rich_time_min;

   /*
    * Timers.
    */
   static tbox::Pointer<tbox::Timer> t_init;
   static tbox::Pointer<tbox::Timer> t_compute_dt;
   static tbox::Pointer<tbox::Timer> t_compute_fluxes;
   static tbox::Pointer<tbox::Timer> t_conservdiff;
   static tbox::Pointer<tbox::Timer> t_setphysbcs;
   static tbox::Pointer<tbox::Timer> t_taggradient;

};

#endif
