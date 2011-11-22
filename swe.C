//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/examples/swe/swe.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Numerical routines for swe equations SAMRAI example
//

#include "swe.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#ifndef LACKS_SSTREAM
#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif
#else
#ifndef included_strstream
#define included_strstream
#include <strstream>
#endif
#endif

using namespace std;

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>

#include "SAMRAI/hier/BoxArray.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

//integer constants for boundary conditions
#define CHECK_BDRY_DATA (0)
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary dat
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"
#include "SAMRAI/appu/CartesianBoundaryUtilities3.h"


// External definitions for Fortran numerical routines
#include "sweFort.h"

// Number of entries in state vector (d_dim veldepth comps +  depth)
#define NSTATE (d_dim.getValue() + 1)

// Problem Spatial Dimension
#define PDIM (d_dim.getValue())

// Number of Scalar Variables
#define NSCAL           (1)

// Number of ghosts cells used for each variable quantity.
#define CELLG           (4)
#define FACEG           (4)
#define FLUXG           (1)

// defines for initialization of particular test cases
// use same values as in cntrl.f90
#define RAMP            (0)
#define DAMBREAKX       (1)
#define DAMBREAKY       (2)
#define ROSSBY          (3)
#define USER_DEFINED    (4)
#define DAMBREAK2D      (5)
#define NOMOTION        (6)
#define NOMOTIONDRYX    (7)
#define THREEHUMP       (8)
#define THREEHUMP_NOMOTION (9)
#define NOMOTIONDRYY    (10)
#define DAMBREAKX_DRY   (11)
#define DAMBREAKX_BERM  (12)
#define THREEHUMPY      (13)
#define DAMBREAKY_DRY   (14)
#define DAMBREAKY_BERM  (15)
#define SAMPSON         (16)
#define TIDETEST        (17)
#define CONRUN          (18)
#define HENICHE         (19)
#define STEP            (20)
#define SUPERCRIT       (21)
#define ROELVINK        (22)
#define ROELVINKY       (23)
#define SLOSH_INLET     (24)
#define TRENCH          (25)

// defines for cell tagging routines
#define RICHARDSON_NEWLY_TAGGED (-10)
#define RICHARDSON_ALREADY_TAGGED (-11)
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif


// Version of swe restart file data
#define SWE_VERSION (1)

// define the various timers
tbox::Pointer<tbox::Timer> swe::t_init;
tbox::Pointer<tbox::Timer> swe::t_compute_dt;
tbox::Pointer<tbox::Timer> swe::t_compute_fluxes;
tbox::Pointer<tbox::Timer> swe::t_conservdiff;
tbox::Pointer<tbox::Timer> swe::t_setphysbcs;
tbox::Pointer<tbox::Timer> swe::t_taggradient;

//======================================================================
//          contructor for swe class
//   Create variables that define the solution vector (w)
//   Set default values or Restart
//   Get from input (which could potentially override restart vals)
//======================================================================

swe::swe(
	const string& object_name,
	const tbox::Dimension& dim,
	tbox::Pointer<tbox::Database> input_db,
	tbox::Pointer<geom::CartesianGridGeometry > grid_geom):
	algs::HyperbolicPatchStrategy(dim),
	d_dim(dim),
	d_nghosts(d_dim),
	d_fluxghosts(d_dim)
{

   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!grid_geom.isNull());

   //set name (string)
   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   //setup timers
   if ( t_init.isNull() ) {
      t_init = tbox::TimerManager::getManager()->
         getTimer("apps::swe::initializeDataOnPatch()");
      t_compute_dt = tbox::TimerManager::getManager()->
         getTimer("apps::swe::computeStableDtOnPatch()");
      t_compute_fluxes = tbox::TimerManager::getManager()->
         getTimer("apps::swe::computeFluxesOnPatch()");
      t_conservdiff = tbox::TimerManager::getManager()->
         getTimer("apps::swe::conservativeDifferenceOnPatch()");
      t_setphysbcs = tbox::TimerManager::getManager()->
         getTimer("apps::swe::setPhysicalBoundaryConditions()");
      t_taggradient = tbox::TimerManager::getManager()->
         getTimer("apps::swe::tagGradientDetectorCells()");
   }

   d_grid_geometry = grid_geom;
  
   d_use_nonuniform_workload = false;

   //define state vector vars w = (h, uh ,vh)^T  
   //define work vars (flux) of size NSTATE 
   d_depth     = new pdat::CellVariable<double>(d_dim, "depth", 1);
   d_drycell   = new pdat::CellVariable<double>(d_dim, "drycell", 1);
   d_addvol    = new pdat::CellVariable<double>(d_dim, "addvol", 1);
   d_bathy     = new pdat::CellVariable<double>(d_dim, "bathymetry",1);
   d_bedlevel  = new pdat::CellVariable<double>(d_dim, "bedlevel",1);
   d_veldepth  = new pdat::CellVariable<double>(d_dim, "veldepth", d_dim.getValue());
   d_fluxm     = new pdat::FaceVariable<double>(d_dim, "fluxm", NSTATE);
   d_fluxp     = new pdat::FaceVariable<double>(d_dim, "fluxp", NSTATE);
  	d_fluxsed   = new pdat::FaceVariable<double>(d_dim, "fluxsed", 1);
   d_scalar    = new pdat::CellVariable<double>(d_dim, "scalar", NSCAL);

   //default physical parameters
   d_rho      = 1025.;  //nominal density for seawater (kg^1/m^3) 
   d_gravity  = 9.81;   //gravitational acceleration (m^1/s^2)  (used for calculating total head for viz)
 
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif

   //number of ghost cells
   d_nghosts    = hier::IntVector(d_dim, CELLG); //cc state vars
   d_fluxghosts = hier::IntVector(d_dim, FLUXG); //flux vars

   //initial data defaults
//   d_radius = tbox::MathUtilities<double>::getSignalingNaN();
//   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_center, d_dim.getValue()); 
//   d_depth_inside = tbox::MathUtilities<double>::getSignalingNaN();

   /*
   d_number_of_intervals = 0;
   d_front_position.resizeArray(0);
   d_interval_density.resizeArray(0);
   d_interval_veldepth.resizeArray(0);
   d_interval_pressure.resizeArray(0);
   */

   //boundary condition arrays - set to bogus values for error check
	if (d_dim == tbox::Dimension(2)){
		d_master_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
		d_scalar_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
		d_vector_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
		for (int ei = 0; ei < NUM_2D_EDGES; ei++) {
			d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
			d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
			d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
		}

		d_master_bdry_node_conds.resizeArray(NUM_2D_NODES);
		d_scalar_bdry_node_conds.resizeArray(NUM_2D_NODES);
		d_vector_bdry_node_conds.resizeArray(NUM_2D_NODES);
		d_node_bdry_edge.resizeArray(NUM_2D_NODES);

		for (int ni = 0; ni < NUM_2D_NODES; ni++) {
			d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
			d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
			d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
			d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
		}

		//These arrays hold the Dirichlet values for bathy/depth/veldepth
		d_bdry_edge_depth.resizeArray(NUM_2D_EDGES);
		d_bdry_edge_bathy.resizeArray(NUM_2D_EDGES);
		d_bdry_edge_bedlevel.resizeArray(NUM_2D_EDGES);
		d_bdry_edge_veldepth.resizeArray(NUM_2D_EDGES*d_dim.getValue());
		tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_edge_depth);
		tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_edge_bedlevel);
		tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_edge_bathy);
		tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_edge_veldepth);
	} //end setup for PDIM=2
	
	//initialize object from restart if restart T
	bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
	if (is_from_restart) {
		getFromRestart();
	}
	//initialize object from input (can overwrite values)
	getFromInput(input_db, is_from_restart);


   //set test problem type from input/restart
   // => d_data_problem_int
   if (d_data_problem == "RAMP") {
      d_data_problem_int = RAMP;  
   } else if (d_data_problem == "ROSSBY") {
      d_data_problem_int = ROSSBY;
   } else if (d_data_problem == "DAMBREAKX") {
      d_data_problem_int = DAMBREAKX;
   } else if (d_data_problem == "DAMBREAKX_DRY") {
      d_data_problem_int = DAMBREAKX_DRY;
   } else if (d_data_problem == "DAMBREAKY_DRY") {
      d_data_problem_int = DAMBREAKY_DRY;
   } else if (d_data_problem == "DAMBREAKX_BERM") {
	  d_data_problem_int = DAMBREAKX_BERM;
   } else if (d_data_problem == "DAMBREAKY_BERM") {
	  d_data_problem_int = DAMBREAKY_BERM;
   } else if (d_data_problem == "DAMBREAKY") { 
	  d_data_problem_int = DAMBREAKY;  
   } else if (d_data_problem == "USER_DEFINED") {
      d_data_problem_int = USER_DEFINED;  
   } else if (d_data_problem == "DAMBREAK2D") {
	  d_data_problem_int = DAMBREAK2D;
   } else if (d_data_problem == "NOMOTION") {
      d_data_problem_int = NOMOTION;
   } else if (d_data_problem == "NOMOTIONDRYX") {
	  d_data_problem_int = NOMOTIONDRYX;
   } else if (d_data_problem == "NOMOTIONDRYY") {
	  d_data_problem_int = NOMOTIONDRYY;
   } else if (d_data_problem == "THREEHUMP") {
      d_data_problem_int = THREEHUMP;
   } else if (d_data_problem == "THREEHUMPY") {
      d_data_problem_int = THREEHUMPY;
   } else if (d_data_problem == "THREEHUMP_NOMOTION") {
	  d_data_problem_int = THREEHUMP_NOMOTION;
   } else if (d_data_problem == "SAMPSON") {
      d_data_problem_int = SAMPSON;
   } else if (d_data_problem == "TIDETEST") {
      d_data_problem_int = TIDETEST;
   } else if (d_data_problem == "CONRUN") {
      d_data_problem_int = CONRUN;
   } else if (d_data_problem == "HENICHE") {
	  d_data_problem_int = HENICHE;
   } else if (d_data_problem == "STEP") {
	  d_data_problem_int = STEP;
   } else if (d_data_problem == "SUPERCRIT") {
	  d_data_problem_int = SUPERCRIT;
	} else if (d_data_problem == "ROELVINK") {
	  d_data_problem_int = ROELVINK;
	} else if (d_data_problem == "ROELVINKY") {
	  d_data_problem_int = ROELVINKY;
	} else if (d_data_problem == "SLOSH_INLET") {
	  d_data_problem_int = SLOSH_INLET;
	} else if (d_data_problem == "TRENCH") {
	  d_data_problem_int = TRENCH;
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "Unknown d_data_problem string = "
         << d_data_problem << " encountered in constructor in swe.C" << endl);
   }

   // set boundary conditions following initialization/restart
	if (d_dim == tbox::Dimension(2)){
   	for (int i = 0; i < NUM_2D_EDGES; i++) {
	      d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
	      d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

	      if (d_master_bdry_edge_conds[i] == REFLECT_BC) {
	         d_scalar_bdry_edge_conds[i] = FLOW_BC;
	      }
	   }

   	for (int i = 0; i < NUM_2D_NODES; i++) {
	      d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
	      d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

	      if (d_master_bdry_node_conds[i] == XREFLECT_BC) {
	         d_scalar_bdry_node_conds[i] = XFLOW_BC;
	      }
	      if (d_master_bdry_node_conds[i] == YREFLECT_BC) {
	         d_scalar_bdry_node_conds[i] = YFLOW_BC;
	      }

	      if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
	         d_node_bdry_edge[i] =
	            appu::CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
	                                            i, d_master_bdry_node_conds[i]);
	      }
	   }
	} // end setting BC for 2D case

   //bind the control parameters between C++ and fortran 
   c2f_(d_dim.getValue(),NSTATE,NSCAL,CELLG,FLUXG, //FORTRAN
	d_data_problem_int,d_C_manning,d_mindepth,d_fluxorder,d_transverse,
	d_sedmodel,d_sedinit,d_taucrit,d_morphfactor); 
              
}  // <=  end of swe constructor

//======================================================================
//          destructor for swe class
//             -nullify the timers
//======================================================================

swe::~swe() 
{
   t_init = NULL;
   t_compute_dt = NULL;
   t_compute_fluxes = NULL;
   t_conservdiff = NULL;
   t_setphysbcs = NULL;
   t_taggradient = NULL;
} // <= end of swe destructor 

//======================================================================
//     register the swe state variables 
//
//       - use a hyperbolic integrator that manages storage
//       - register the variable type (TIME_DEP , FLUX, etc.)
//       - register the data with the visualization 
//       - specify intrinsic or user-defined functions for restr/prolong
//       - coarsening of depth is done using standard SAMRAI libs
//       - coarsening of veldepth is done using standard SAMRAI libs
//
//     this is a concrete implementation of a function from
//     the algs::HyperbolicPatchStrategy abstract base class
//======================================================================

void swe::registerModelVariables(algs::HyperbolicLevelIntegrator* integrator)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(integrator != (algs::HyperbolicLevelIntegrator *)NULL);
   TBOX_ASSERT(CELLG == FACEG);
#endif

   // register state variables to work in AMR scheme
   integrator->registerVariable(
		d_depth ,d_nghosts, 
		algs::HyperbolicLevelIntegrator::TIME_DEP,
		d_grid_geometry,
		"CONSERVATIVE_COARSEN",
		"CONSERVATIVE_LINEAR_REFINE");

   integrator->registerVariable(
		d_bathy ,d_nghosts, 
		algs::HyperbolicLevelIntegrator::TIME_DEP,
		d_grid_geometry,
		"CONSERVATIVE_COARSEN",
		"CONSERVATIVE_LINEAR_REFINE");
		
	integrator->registerVariable(
		d_bedlevel ,d_nghosts, 
		algs::HyperbolicLevelIntegrator::TIME_DEP,
		d_grid_geometry,
		"CONSERVATIVE_COARSEN",
		"CONSERVATIVE_LINEAR_REFINE");

   integrator->registerVariable(
		d_drycell ,d_nghosts, 
		algs::HyperbolicLevelIntegrator::TIME_DEP,
		d_grid_geometry,
		"NO_COARSEN",
		"NO_REFINE");

   integrator->registerVariable(
		d_addvol ,d_nghosts, 
		algs::HyperbolicLevelIntegrator::TIME_DEP,
		d_grid_geometry,
		"CONSERVATIVE_COARSEN",
		"CONSERVATIVE_LINEAR_REFINE");

   integrator->registerVariable(
		d_veldepth, d_nghosts, 
		algs::HyperbolicLevelIntegrator::TIME_DEP,
		d_grid_geometry,
		"CONSERVATIVE_COARSEN",
		"CONSERVATIVE_LINEAR_REFINE");

   integrator->registerVariable(
		d_scalar, d_nghosts, 
		algs::HyperbolicLevelIntegrator::TIME_DEP,
		d_grid_geometry,
		"CONSERVATIVE_COARSEN",
		"CONSERVATIVE_LINEAR_REFINE");
                               
   integrator->registerVariable(
		d_fluxm, d_fluxghosts, 
		algs::HyperbolicLevelIntegrator::FLUX,
		d_grid_geometry,
		"CONSERVATIVE_COARSEN",
		"NO_REFINE");

   integrator->registerVariable(
		d_fluxp, d_fluxghosts, 
		algs::HyperbolicLevelIntegrator::FLUX,
		d_grid_geometry,
		"CONSERVATIVE_COARSEN",
		"NO_REFINE");
		
	integrator->registerVariable(
		d_fluxsed, d_fluxghosts, 
		algs::HyperbolicLevelIntegrator::FLUX,
		d_grid_geometry,
		"CONSERVATIVE_COARSEN",
		"NO_REFINE");


   //the variable database is used for setting up output (visualization)
   //is then discarded
   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

   d_plot_context = integrator->getPlotContext();

   //visit setup
#ifdef HAVE_HDF5   
   if (!(d_visit_writer.isNull())) {
	
		//prognostic vars
      d_visit_writer->registerPlotQuantity(
			"depth","SCALAR",vardb->mapVariableAndContextToIndex(d_depth, d_plot_context));
			
      d_visit_writer->registerPlotQuantity(
			"bathy","SCALAR",vardb->mapVariableAndContextToIndex(d_bathy, d_plot_context));
	  
		d_visit_writer->registerPlotQuantity(
			"bedlevel","SCALAR",vardb->mapVariableAndContextToIndex(d_bedlevel, d_plot_context));
			
		d_visit_writer->registerPlotQuantity(
			"drycell","SCALAR", vardb->mapVariableAndContextToIndex( d_drycell, d_plot_context));
			    
      d_visit_writer->registerPlotQuantity(
			"VD", "VECTOR", vardb->mapVariableAndContextToIndex( d_veldepth, d_plot_context));
			
		//diagnostic vars
      d_visit_writer->registerDerivedPlotQuantity("velocity", "VECTOR", this);
			
      d_visit_writer->registerDerivedPlotQuantity("Total Energy", "SCALAR", this);
			
      d_visit_writer->registerDerivedPlotQuantity("u", "SCALAR", this);				
									
		d_visit_writer->registerDerivedPlotQuantity("v", "SCALAR", this);
		
		d_visit_writer->registerDerivedPlotQuantity("zeta", "SCALAR", this);
						
		d_visit_writer->registerDerivedPlotQuantity("wetdry", "SCALAR", this);
													  																							  
   }

   if (d_visit_writer.isNull()) {
      TBOX_WARNING(d_object_name << ": registerModelVariables()\n"
		   << "A Visit data writer was\n"
		   << "registered.  Consequently, no plot data will\n"
		   << "be written." << endl);
   }
#endif
	 
} // <= end registering the variables and setting up visualization

//======================================================================
//     setup nonuniform load balancing (may not be used)
//======================================================================

void swe::setupLoadBalancer(
	algs::HyperbolicLevelIntegrator* integrator,
	mesh::GriddingAlgorithm* gridding_algorithm)
	{
   (void)integrator;

   const hier::IntVector& zero_vec = hier::IntVector::getZero(d_dim);

   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

   if (d_use_nonuniform_workload && gridding_algorithm) {
      tbox::Pointer<mesh::TreeLoadBalancer> load_balancer =
         gridding_algorithm->getLoadBalanceStrategy();

      if (!load_balancer.isNull()) {
         d_workload_variable = new pdat::CellVariable<double>(
               d_dim,
               "workload_variable",
               1);
         d_workload_data_id =
            vardb->registerVariableAndContext(d_workload_variable,
               vardb->getContext("WORKLOAD"),
               zero_vec);
         load_balancer->setWorkloadPatchDataIndex(d_workload_data_id);
         vardb->registerPatchDataForRestart(d_workload_data_id);
      } else {
         TBOX_WARNING(
            d_object_name << ": "
            <<
            "  Unknown load balancer used in gridding algorithm."
            <<
            "  Ignoring request for nonuniform load balancing." << endl);
         d_use_nonuniform_workload = false;
      }
   } else {
      d_use_nonuniform_workload = false;
   }

} // <= end setup of non-uniform load balance 

//======================================================================
//    set initial state on a patch (no restart) 
//       note: during sim, new patches get data through interpolation
//             ----> this routine is initial model state only
//
//     this is a concrete implementation of a function from
//     the algs::HyperbolicPatchStrategy abstract base class
//======================================================================

void swe::initializeDataOnPatch(
	hier::Patch& patch,
   const double data_time,
	const bool initial_time)
{
   (void)data_time;

   t_init->start();

   if (initial_time) {

      const tbox::Pointer<geom::CartesianPatchGeometry > pgeom = patch.getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();
      const double* xhi = pgeom->getXUpper();

      //pointer to patch data h and uh    
      tbox::Pointer< pdat::CellData<double> > depth  = 
          patch.getPatchData(d_depth, getDataContext());
      tbox::Pointer< pdat::CellData<double> > veldepth = 
          patch.getPatchData(d_veldepth, getDataContext());
      tbox::Pointer< pdat::CellData<double> > bathy = 
          patch.getPatchData(d_bathy, getDataContext());
      tbox::Pointer< pdat::CellData<double> > addvol = 
          patch.getPatchData(d_addvol, getDataContext());
      tbox::Pointer< pdat::CellData<double> > drycell = 
          patch.getPatchData(d_drycell, getDataContext());
      tbox::Pointer< pdat::CellData<double> > bedlevel = 
          patch.getPatchData(d_bedlevel, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!depth.isNull());
      TBOX_ASSERT(!veldepth.isNull());
      TBOX_ASSERT(!bedlevel.isNull());
#endif
      const hier::IntVector& ghost_cells = veldepth->getGhostCellWidth();
	
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(veldepth->getGhostCellWidth() == ghost_cells);
#endif

      const hier::Index ifirst=patch.getBox().lower();
      const hier::Index ilast =patch.getBox().upper();

      //initialize the volume add to zero
     //addvol->fillAll(0);

      //set bathymetry 
      setbathy_(d_data_problem_int, //FORTRAN
                   dx,xlo,xhi,
                   ifirst(0),ilast(0),
                   ifirst(1),ilast(1),
                   ghost_cells(0),
                   ghost_cells(1),
                   bathy->getPointer());

      //initialize based on choice of problem
      initflow_(d_data_problem_int, //FORTRAN
                   dx,xlo,xhi,
                   ifirst(0),ilast(0),
                   ifirst(1),ilast(1),
                   ghost_cells(0),
                   ghost_cells(1),
                   depth->getPointer(),
                   veldepth->getPointer(),
                   bathy->getPointer(),
                   bedlevel->getPointer());

   } //if initial_time

   if (d_use_nonuniform_workload) {
      if (!patch.checkAllocated(d_workload_data_id)) {
         patch.allocatePatchData(d_workload_data_id);
      }
      tbox::Pointer< pdat::CellData<double> > workload_data =
         patch.getPatchData(d_workload_data_id);
      workload_data->fillAll(1.0);
   }

   t_init->stop();


} // <= end set initial conditions 

//======================================================================
//    compute stable time step on a patch 
//      returns stable DT (double)
//
//     this is a concrete implementation of a function from
//     the algs::HyperbolicPatchStrategy abstract base class
//======================================================================

double swe::computeStableDtOnPatch(
   hier::Patch& patch,
   const bool initial_time, 
   const double dt_time) 
{
   (void) initial_time;
   (void) dt_time;

   //start timer
   t_compute_dt->start();

   //grab delta_x from the patch geometry
   const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();
 
   //grab computational dimensions of the patch
   const hier::Index ifirst=patch.getBox().lower();
   const hier::Index ilast =patch.getBox().upper();

   //grab pointer to state vector [w] 
   const tbox::Pointer< pdat::CellData<double> > depth    = 
      patch.getPatchData(d_depth, getDataContext());
   const tbox::Pointer< pdat::CellData<double> > veldepth = 
      patch.getPatchData(d_veldepth, getDataContext());
   const tbox::Pointer< pdat::CellData<double> > bathy    = 
      patch.getPatchData(d_bathy, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!depth.isNull());
   TBOX_ASSERT(!veldepth.isNull());
#endif 
   const hier::IntVector& ghost_cells = depth->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(veldepth->getGhostCellWidth() == ghost_cells);
#endif 

   //call routine to calculate min deltaT
   //not sure why, but ghost cells are included in the DT calculation
   double stabdt = 0.;
   calcdt_(dx,  //FORTRAN
             ifirst(0),ilast(0),
   			 ifirst(1),ilast(1),
   			 ghost_cells(0),
   			 ghost_cells(1),
             depth->getPointer(),
             veldepth->getPointer(),
             bathy->getPointer(),
             stabdt);
   //cout << "time step" << stabdt; 
   //stop the timer
   t_compute_dt->stop();
   return stabdt;

} // <= end calculate stable time step on a patch

//======================================================================
//    compute fluxes at each cell face on a patch 
//
//     this is a concrete implementation of a function from
//     the algs::HyperbolicPatchStrategy abstract base class
//======================================================================


void swe::computeFluxesOnPatch(
	hier::Patch& patch,
	const double time, 
	const double dt)
{
   (void) time;

   //start timing of flux computation
   t_compute_fluxes->start();


#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif

   //set patch_geometry + dx
   const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();
   const double* xhi = patch_geom->getXUpper();

   //set computational limits
   hier::Box pbox = patch.getBox();
   const hier::Index ifirst= pbox.lower();
   const hier::Index ilast = pbox.upper();

   //get pointers to H, vH, flux
   tbox::Pointer< pdat::CellData<double> > depth    = 
      patch.getPatchData(d_depth, getDataContext());
   tbox::Pointer< pdat::CellData<double> > bathy    = 
      patch.getPatchData(d_bathy, getDataContext());
   tbox::Pointer< pdat::CellData<double> > bedlevel = 
      patch.getPatchData(d_bedlevel, getDataContext());
   tbox::Pointer< pdat::CellData<double> > veldepth = 
      patch.getPatchData(d_veldepth, getDataContext());
   tbox::Pointer< pdat::FaceData<double> > fluxm    = 
      patch.getPatchData(d_fluxm, getDataContext());
   tbox::Pointer< pdat::FaceData<double> > fluxp    = 
       patch.getPatchData(d_fluxp, getDataContext());
   tbox::Pointer< pdat::FaceData<double> > fluxsed    = 
      patch.getPatchData(d_fluxsed, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!bathy.isNull());
	TBOX_ASSERT(!depth.isNull());
	TBOX_ASSERT(!bedlevel.isNull());
   TBOX_ASSERT(!veldepth.isNull());
   TBOX_ASSERT(!fluxm.isNull());
	TBOX_ASSERT(!fluxp.isNull());
	TBOX_ASSERT(!fluxsed.isNull());
   TBOX_ASSERT(depth->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(bedlevel->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(veldepth->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(bathy->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(fluxm->getGhostCellWidth() == d_fluxghosts);
   TBOX_ASSERT(fluxp->getGhostCellWidth() == d_fluxghosts);
   TBOX_ASSERT(fluxsed->getGhostCellWidth() == d_fluxghosts);
#endif

   
   // D.L. George scheme 
   // Inputs  - h,vh,b 
   // Outputs - flux_x , flux_y 
   george_flux_(d_data_problem_int,dt,dx, //FORTRAN
              ifirst(0),ilast(0),ifirst(1),ilast(1),
              depth->getPointer(),
              veldepth->getPointer(),
              bathy->getPointer(),
              fluxm->getPointer(0),
   			  fluxm->getPointer(1),
   			  fluxp->getPointer(0),
   		     fluxp->getPointer(1));

   if (d_sedmodel > 0 && time > d_sedinit) {
		flux_sed_(d_data_problem_int,dt,dx,xlo,xhi, //FORTRAN
			ifirst(0),ilast(0),ifirst(1),ilast(1),
			depth->getPointer(),
			veldepth->getPointer(),
			bedlevel->getPointer(),
			fluxsed->getPointer(0),
			fluxsed->getPointer(1));
	 }

   //stop flux computation timer
   t_compute_fluxes->stop();
}
//======================================================================
//          conservativeDifferenceOnPatch
//
//    update solution variables to new time level using flux difference
//
//    this is a concrete implementation of a function from
//    the algs::HyperbolicPatchStrategy abstract base class
//
//    note the fluxes are already integrated over deltaT 
//    so no time step is used in the cons. difference
//    but given that it is passed in the header, probably could
//    be used here if fluxes were not integrated over deltaT 
//
//    this is a concrete implementation of a function from
//    the algs::HyperbolicPatchStrategy abstract base class
//
//     in this routine the fortran90 function consdiff is called   
//     to perform the flux difference.  This function first converts 
//     to conservative form, summs the fluxes (actually the difference)
//     then converts back to primitive
//======================================================================

void swe::conservativeDifferenceOnPatch(
	hier::Patch& patch,
	const double time,
	const double dt,
	bool at_syncronization)
{
   (void) time;
   (void) dt;
   (void) at_syncronization;

   //start the timer
   t_conservdiff->start();

   //get patch mesh spacing (delta_x,delta_y)
   const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();
   const double* xhi = patch_geom->getXUpper();

   //get patch array limits (i1,i2,j1,j2)
   hier::Box pbox = patch.getBox();
   const hier::Index ifirst= pbox.lower();
   const hier::Index ilast = pbox.upper();
 
   //get pointers to model fields
   tbox::Pointer< pdat::CellData<double> > depth       = 
      patch.getPatchData(d_depth, getDataContext());
   tbox::Pointer< pdat::CellData<double> > bathy       = 
      patch.getPatchData(d_bathy, getDataContext());
   tbox::Pointer< pdat::CellData<double> > bedlevel    = 
      patch.getPatchData(d_bedlevel, getDataContext());
   tbox::Pointer< pdat::CellData<double> > veldepth    = 
      patch.getPatchData(d_veldepth, getDataContext());
   tbox::Pointer< pdat::FaceData<double> > fluxm       = 
      patch.getPatchData(d_fluxm, getDataContext());
   tbox::Pointer< pdat::FaceData<double> > fluxp       = 
      patch.getPatchData(d_fluxp, getDataContext());
   tbox::Pointer< pdat::FaceData<double> > fluxsed       = 
      patch.getPatchData(d_fluxsed, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!depth.isNull());
   TBOX_ASSERT(!veldepth.isNull());
   TBOX_ASSERT(!bedlevel.isNull());
   TBOX_ASSERT(!fluxm.isNull());
	TBOX_ASSERT(!fluxp.isNull());
	TBOX_ASSERT(!fluxsed.isNull());
   TBOX_ASSERT(depth->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(veldepth->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(fluxm->getGhostCellWidth() == d_fluxghosts);
   TBOX_ASSERT(fluxp->getGhostCellWidth() == d_fluxghosts);
   TBOX_ASSERT(fluxsed->getGhostCellWidth() == d_fluxghosts);
#endif

   if (d_dim == tbox::Dimension(2)) {
			
			
		// update state variables using summation of fluxes
		consdiff_(dx,ifirst(0),ilast(0),ifirst(1),ilast(1), //FORTRAN
			fluxm->getPointer(0),
			fluxm->getPointer(1),
			fluxp->getPointer(0),
			fluxp->getPointer(1),
			depth->getPointer(),
			veldepth->getPointer(),
			bathy->getPointer());

		//update bed thickness if sediment model is active and bathymetry if morpho active
		if(d_sedmodel > 0 && time > d_sedinit){
			consdiff_sed_(dx,ifirst(0),ilast(0),ifirst(1),ilast(1), //FORTRAN
				fluxsed->getPointer(0),
				fluxsed->getPointer(1),
				bedlevel->getPointer(),
				bathy->getPointer());
		}
		
		// friction term (d_frictype==0 => no friction)
		if (d_frictype == 1) { //implicit Linear friction 
			linfriction_(dx,dt,ifirst(0),ilast(0),ifirst(1),ilast(1), //FORTRAN
				depth->getPointer(),
				veldepth->getPointer(), 
				bathy->getPointer());
		} else if(d_frictype == 2) { //implicit Manning friction with coefficinet n=C_manning
			friction_(dx,dt,ifirst(0),ilast(0),ifirst(1),ilast(1), //FORTRAN
				depth->getPointer(),
				veldepth->getPointer(), 
				bathy->getPointer());
		}

		// dry check - make sure cells have non-negative depth 
		drycheck_(dx,dt,ifirst(0),ilast(0),ifirst(1),ilast(1), //FORTRAN
			depth->getPointer(),
			veldepth->getPointer()); 

		// junk probe 
		junkprobe_(d_data_problem_int,dx,xlo,xhi,time, //FORTRAN
			ifirst(0),ilast(0),ifirst(1),ilast(1),
			depth->getPointer(),
			veldepth->getPointer(),
			bathy->getPointer());
	}
	
   //stop timer
   t_conservdiff->stop();

} // <-- end conservative differencing

//=======================================================================
//       boundaryReset      
//
//  reset physical boundary conditions for reflective 
//=======================================================================

void swe::boundaryReset(
	hier::Patch& patch,
	pdat::FaceData<double>& traced_left,  
	pdat::FaceData<double>& traced_right) const
{
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();
   int i, idir;
   bool bdry_cell = true;

   const tbox::Pointer<geom::CartesianPatchGeometry> patch_geom =
      patch.getPatchGeometry();
   hier::BoxArray domain_boxes(d_dim);
   d_grid_geometry->computePhysicalDomain(domain_boxes, patch_geom->getRatio());
   int num_domain_boxes = domain_boxes.getNumberOfBoxes();
   const double* dx = patch_geom->getDx();
   const double* xpatchhi = patch_geom->getXUpper();
   const double* xdomainhi = d_grid_geometry->getXUpper();
   const double* xpatchlo = patch_geom->getXLower();
   const double* xdomainlo = d_grid_geometry->getXLower();

   pdat::CellIndex icell(ifirst);
   hier::BoxArray bdrybox(d_dim, 2 * d_dim.getValue());
   hier::Index ibfirst = ifirst;
   hier::Index iblast = ilast;
   int bdry_case = 0;

   for (idir = 0; idir < d_dim.getValue(); idir++) {
      ibfirst(idir) = ifirst(idir) - 1;
      iblast(idir) = ifirst(idir) - 1;
      bdrybox[2 * idir] = hier::Box(ibfirst, iblast);

      ibfirst(idir) = ilast(idir) + 1;
      iblast(idir) = ilast(idir) + 1;
      bdrybox[2 * idir + 1] = hier::Box(ibfirst, iblast);
   }

   for (idir = 0; idir < d_dim.getValue(); idir++) {
      int bside = 2 * idir;
      if (d_dim == tbox::Dimension(2)) {
         bdry_case = d_master_bdry_edge_conds[bside];
      }
      if (bdry_case == REFLECT_BC) {
         for (pdat::CellIterator ic(bdrybox[bside]); ic; ic++) {
            for (i = 0; i < num_domain_boxes; i++) {
               if (domain_boxes[i].contains(ic()))
                  bdry_cell = false;
            }
            if (bdry_cell) {
               pdat::FaceIndex sidein = pdat::FaceIndex(ic(), idir, 1);
               (traced_left)(sidein, 0) = (traced_right)(sidein, 0);
            }
         }
      }

      int bnode = 2 * idir + 1;
     	 if (d_dim == tbox::Dimension(2)) {
	         bdry_case = d_master_bdry_edge_conds[bnode];
	      }
      
// BEGIN SIMPLE-MINDED FIX FOR STEP PROBLEM
      if ((d_data_problem == "STEP") && (bnode == 1) &&
          (tbox::MathUtilities<double>::Abs(xpatchhi[0] - xdomainhi[0]) < dx[0])) {
         bdry_case = FLOW_BC;
      }
// // BEGIN SIMPLE-MINDED FIX FOR ROELVINK PROBLEM
      if ((d_data_problem == "ROELVINK") && (bnode == 1) &&
           (tbox::MathUtilities<double>::Abs(xpatchlo[0] - xdomainlo[0]) < dx[0])) {
          bdry_case = DIRICHLET_BC;
       }
// END SIMPLE-MINDED FIX FOR ROELVINK PROBLEM
// // BEGIN SIMPLE-MINDED FIX FOR ROELVINKY PROBLEM
      if ((d_data_problem == "ROELVINKY") && (bnode == 3) &&
           (tbox::MathUtilities<double>::Abs(xpatchlo[1] - xdomainlo[1]) < dx[1])) {
          bdry_case = DIRICHLET_BC;
       }
// END SIMPLE-MINDED FIX FOR ROELVINKY PROBLEM
      if (bdry_case == REFLECT_BC) {
         for (pdat::CellIterator ic(bdrybox[bnode]); ic; ic++) {
            for (i = 0; i < num_domain_boxes; i++) {
               if (domain_boxes[i].contains(ic()))
                  bdry_cell = false;
            }
            if (bdry_cell) {
               pdat::FaceIndex sidein = pdat::FaceIndex(ic(), idir, 0);
               (traced_right)(sidein, 0) = (traced_left)(sidein, 0);
            }
         }
      }
   }
}

//======================================================================
//    setPhysicalBoundaryConditions
//
//  set data in the ghost cells corresponding to physical boundary 
//  conditions. 
//======================================================================

void swe::setPhysicalBoundaryConditions(
   hier::Patch& patch, 
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   (void) fill_time;

   //start time
   t_setphysbcs->start();

   //get pointer to state vars
   tbox::Pointer< pdat::CellData<double> > depth  =
      patch.getPatchData(d_depth, getDataContext());
   tbox::Pointer< pdat::CellData<double> > veldepth =
      patch.getPatchData(d_veldepth, getDataContext());
   tbox::Pointer< pdat::CellData<double> > bathy = 
      patch.getPatchData(d_bathy, getDataContext());
   tbox::Pointer< pdat::CellData<double> > bedlevel = 
      patch.getPatchData(d_bedlevel, getDataContext());


#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!depth.isNull());
   TBOX_ASSERT(!veldepth.isNull());
   TBOX_ASSERT(!bedlevel.isNull());
#endif
   hier::IntVector ghost_cells = depth->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(veldepth->getGhostCellWidth() == ghost_cells);
#endif


   //set boundary conditions for cells corresponding to patch edges
   tbox::Array<int> tmp_edge_scalar_bcond(NUM_2D_EDGES);
   tbox::Array<int> tmp_edge_veldepth_bcond(NUM_2D_EDGES);
   tbox::Array<int> tmp_edge_vector_bcond(NUM_2D_EDGES);
   for (int i = 0; i < NUM_2D_EDGES; i++) {
      tmp_edge_scalar_bcond[i] = d_scalar_bdry_edge_conds[i];
      tmp_edge_veldepth_bcond[i] = d_scalar_bdry_edge_conds[i];
      tmp_edge_vector_bcond[i] = d_vector_bdry_edge_conds[i];
   }

	//explicitly set time-dependent elevation for certain test cases
  if (d_data_problem == "RAMP") {
    for (int i = 0; i < NUM_2D_EDGES; i++) {
      if(d_scalar_bdry_edge_conds[i] == DIRICHLET_BC){   
			d_bdry_edge_depth[i] = 5.+cos(fill_time*2*3.1415/43920.);
			//cout << "bc time " << fill_time;
  		}
    }
	} else if(d_data_problem == "CONRUN"){
		double eta0_cr = .32;
		double Hw_cr = .032;
		double cs_cr = sqrt(9.81*eta0_cr)*(1 + Hw_cr/(2*eta0_cr));
		double ls_cr = eta0_cr*sqrt( (4*cs_cr*eta0_cr)/(3*Hw_cr*sqrt(9.81*eta0_cr)));
		double T_cr  = 2.45;
		double u, eta;
		for (int i = 0; i < NUM_2D_EDGES; i++) {
			if(d_scalar_bdry_edge_conds[i] == DIRICHLET_BC){   
			//equation 30, valiani and begnudelli, JHE,2006	
				eta = Hw_cr/pow( cosh(cs_cr*(fill_time-T_cr)/ls_cr),2.0);
			//equation 3 , Choi et al. coastal engineering, v54, 2007, 
				u = .75*((d_bdry_edge_depth[i]-eta0_cr)/eta0_cr)*sqrt(9.81*(eta0_cr+Hw_cr));
				d_bdry_edge_depth[i] = eta0_cr + eta;
				d_bdry_edge_veldepth[i*PDIM+0] = u*d_bdry_edge_depth[i];
			}
		}
	} else if(d_data_problem == "HENICHE"){
		double H0_hen = 1.0;
		double eta0_hen = 0.75;
		double T_hen = 3600;
		double h; //,u;
		for (int i = 0; i < NUM_2D_EDGES; i++) {
			if(d_scalar_bdry_edge_conds[i] == DIRICHLET_BC){   
			//equation 33, Brufau et al, IJNMF v39
				h = H0_hen + eta0_hen*cos(2*3.14159*fill_time/T_hen);
			//equation  
			//u = eta0_hen*sqrt(9.81/H0_hen)*cos(2*3.14159*fill_time/T_hen);
				d_bdry_edge_depth[i] = h;
			//d_bdry_edge_veldepth[i*PDIM+0] = u*h;
			}
		}
		// set open boundary forcing for ROELVINK test where open boundary is at X=0 
		// Note we are forcing the free surface (actually the depth) but do not specify the velocity
		// thus we are resetting the scalar_bcond to DIRICHLET (which then uses the values we provide)
		// and we are resetting the vector_bcond to FLOW (which seems to automatically use zeroth-order 
		// extrapolation for ud/vd)
		// note we do not set values for d_bdry_edge_veldepth
   } else if(d_data_problem == "ROELVINK"){
		double H0_roel = 2;
		double eta0_roel = 1.0;
		double T_roel = 12*3600;
		double h; //,u;
		const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom = patch.getPatchGeometry();
		const double* dx = patch_geom->getDx();
		const double* xpatchlow = patch_geom->getXLower();
		const double* xdomainlow = d_grid_geometry->getXLower();
		if (tbox::MathUtilities<double>::Abs(xpatchlow[0]-xdomainlow[0]) < dx[0]) {
			tmp_edge_scalar_bcond[XLO] = DIRICHLET_BC;
		 	tmp_edge_vector_bcond[XLO] = FLOW_BC;
		   tmp_edge_veldepth_bcond[XLO] = FLOW_BC;
		}
		for (int i = 0; i < NUM_2D_EDGES; i++) {
										if(tmp_edge_scalar_bcond[i] == DIRICHLET_BC){
			//if(d_scalar_bdry_edge_conds[i] == DIRICHLET_BC){
										//equation 33, Brufau et al, IJNMF v39
				h = H0_roel + eta0_roel*cos(2*3.14159*fill_time/T_roel);
				//cout << "setting h on edge " << h << "  " << i << '\n';
										//equation  
										//u = eta0_hen*sqrt(9.81/H0_hen)*cos(2*3.14159*fill_time/T_hen);
				d_bdry_edge_depth[i] = h;
				d_bdry_edge_bathy[i] = -2;  //MUST SET DIRICHLET BATHYMETRY
				//d_bdry_edge_bedlevel[i] = 99.0;
				//								d_bdry_edge_veldepth[i*PDIM+0] = 0.; //u*h;
				//								d_bdry_edge_veldepth[i*PDIM+1] = 0.; //u*h;
			}
		}
	// set open boundary forcing for ROELVINK test where open boundary is at Y=0. 
	// Note we are forcing the free surface (actually the depth) but do not specify the velocity
	// thus we are resetting the scalar_bcond to DIRICHLET (which then uses the values we provide)
	// and we are resetting the vector_bcond to FLOW (which seems to automatically use zeroth-order 
	// extrapolation for ud/vd)
	// note we do not set values for d_bdry_edge_veldepth 
	} else if(d_data_problem == "ROELVINKY"){
		double H0_roel = 10;
		double eta0_roel = 1.0;
		double T_roel = 12*3600;
		double h; //,u;
		const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom = patch.getPatchGeometry();
		const double* dx = patch_geom->getDx();
		const double* xpatchlow = patch_geom->getXLower();
		const double* xdomainlow = d_grid_geometry->getXLower();
		if (tbox::MathUtilities<double>::Abs(xpatchlow[1]-xdomainlow[1]) < dx[1]) {
			tmp_edge_scalar_bcond[YLO] = DIRICHLET_BC;
			tmp_edge_vector_bcond[YLO] = FLOW_BC;
		}
		for (int i = 0; i < NUM_2D_EDGES; i++) {
										if(tmp_edge_scalar_bcond[i] == DIRICHLET_BC){
			//if(d_scalar_bdry_edge_conds[i] == DIRICHLET_BC){
										//equation 33, Brufau et al, IJNMF v39
				h = H0_roel + eta0_roel*cos(2*3.14159*fill_time/T_roel);
				//cout << "setting h on edge " << h << "  " << i << '\n';
										//equation  
										//u = eta0_hen*sqrt(9.81/H0_hen)*cos(2*3.14159*fill_time/T_hen);
				d_bdry_edge_depth[i] = h;
				d_bdry_edge_bathy[i] = -10;  //MUST SET DIRICHLET BATHYMETRY
												//d_bdry_edge_veldepth[i*PDIM+0] = 0.; //u*h;
												//d_bdry_edge_veldepth[i*PDIM+1] = 0.; //u*h;
			}
		}
		
		
	} //end if d_data_problem
	
	//special fix for step, set flow conditions only if it is truly the exit 
	if (d_data_problem == "STEP") {

      const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom =
         patch.getPatchGeometry();
      const double* dx = patch_geom->getDx();
      const double* xpatchhi = patch_geom->getXUpper();
      const double* xdomainhi = d_grid_geometry->getXUpper();

      if (tbox::MathUtilities<double>::Abs(xpatchhi[0]-xdomainhi[0]) < dx[0]) {
         tmp_edge_scalar_bcond[XHI] = FLOW_BC;
         tmp_edge_vector_bcond[XHI] = FLOW_BC;
      }

   }

 //        //special fix for step, set dirichlet conditions only if it is truly the open boundary
 //         if (d_data_problem == "ROELVINK") {
 // 
 //       const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom =
 //          patch.getPatchGeometry();
 //       const double* dx = patch_geom->getDx();
 //       const double* xpatchlow = patch_geom->getXLower();
 //       const double* xdomainlow = d_grid_geometry->getXLower();
 // //gwc
 //       if (tbox::MathUtilities<double>::Abs(xpatchlow[0]-xdomainlow[0]) > dx[0]) {
 //          tmp_edge_scalar_bcond[XLO] = DIRICHLET_BC;
 //          tmp_edge_vector_bcond[XLO] = DIRICHLET_BC;
 //       }
 //    }


   


   appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("depth", depth,
                           patch,
                           ghost_width_to_fill,
                           tmp_edge_scalar_bcond,
                           d_bdry_edge_depth);
 
   //without this, threehump_nomotion fails on the boundaries
   appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("bathy", bathy,
   	                   patch,
   					       ghost_width_to_fill,
   					       tmp_edge_scalar_bcond,
   					       d_bdry_edge_bathy);

	appu::CartesianBoundaryUtilities2::
		fillEdgeBoundaryData("bedlevel", bedlevel,
			patch,
			ghost_width_to_fill,
			tmp_edge_veldepth_bcond,
			d_bdry_edge_bedlevel);
					
   appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("veldepth", veldepth,
                           patch,
                           ghost_width_to_fill,
                           tmp_edge_vector_bcond,
                           d_bdry_edge_veldepth);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(EDGE2D_BDRY_TYPE, patch, ghost_width_to_fill,
                     tmp_edge_scalar_bcond, tmp_edge_vector_bcond);
#endif
#endif

   
   //set boundary conditions for cells corresponding to patch nodes
   appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("depth", depth,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_edge_depth);

   appu::CartesianBoundaryUtilities2::
	  fillNodeBoundaryData("depth", bathy,
					       patch,
					       ghost_width_to_fill,
					       d_scalar_bdry_node_conds,
					       d_bdry_edge_bathy);
					
	appu::CartesianBoundaryUtilities2::
		  fillNodeBoundaryData("bedlevel", bedlevel,
						       patch,
						       ghost_width_to_fill,
						       d_scalar_bdry_node_conds,
						       d_bdry_edge_bedlevel);
					
   appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("veldepth", veldepth,
                           patch,
                           ghost_width_to_fill,
                           d_vector_bdry_node_conds,
                           d_bdry_edge_veldepth);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(NODE2D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_node_conds, d_vector_bdry_node_conds);
#endif
#endif
			

   // end timer
   t_setphysbcs->stop();
}

//======================================================================
//    tagGradientDetectorCells
//
//    Tag cells for refinement using gradient detector using criteria
//       defined in the input
//   
//    Operates: patch by patch
//
//    this is a concrete implementation of a function from
//    the algs::HyperbolicPatchStrategy<PDIM> abstract base class
//======================================================================

void swe::tagGradientDetectorCells(
	hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_indx,
   const bool uses_richardson_extrapolation_too)
{
   (void) initial_error;

   //start timer 
   t_taggradient->start();

   const int error_level_number = patch.getPatchLevelNumber();

   //get patch geometry, dx/dy, and the bounding box
   const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom = patch.getPatchGeometry();
	const double* dx  = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();
   const double* xhi = patch_geom->getXUpper();
 
   tbox::Pointer< pdat::CellData<int> > tags        = patch.getPatchData(tag_indx);

   hier::Box pbox = patch.getBox(); 
   //hier::Box pboxm1 = pbox.grow(pbox,-1);
   hier::BoxArray domain_boxes(d_dim);
   d_grid_geometry->computePhysicalDomain(domain_boxes, patch_geom->getRatio());
   /*
    * Construct domain bounding box
    */
   hier::Box domain(d_dim);
   for( int i=0; i < domain_boxes.getNumberOfBoxes(); i++ ) {
     domain += domain_boxes[i];
   }

   const hier::Index domfirst=domain.lower();
   const hier::Index domlast =domain.upper();
   const hier::Index ifirst=patch.getBox().lower();
   const hier::Index ilast =patch.getBox().upper();

//cout << "setting tags on patch  " << xlo[0] <<"\t"<< xlo[1] << "\t" << xhi[0] << "\t" << xhi[1] << "\t" << error_level_number << "\n";
 
   //allocate temporary tags and set to not tagged
   tbox::Pointer<pdat::CellData<int> > temp_tags(new pdat::CellData<int>(
                                                    pbox,
                                                    1,
                                                    d_nghosts));
   temp_tags->fillAll(FALSE);


   /*
    * Possible tagging criteria includes 
    *    DEPTH_GRADIENT 
    *    BATHY_GRADIENT
	 *    U_GRADIENT
    * The criteria is specified over a time interval.
    *
    * Loop over criteria provided and check to make sure we are in the
    * specified time interval.  If so, apply appropriate tagging for
    * the level.
    */
   for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ncrit++) {

      string ref = d_refinement_criteria[ncrit];
      tbox::Pointer< pdat::CellData<double > > var;     
      int size = 0;
      double tol = 0.;
      bool time_allowed = false;
      //tag based on DEPTH_GRADIENT
      if (ref == "DEPTH_GRADIENT") {
         var = patch.getPatchData(d_depth, getDataContext());
         size = d_depth_grad_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_depth_grad_tol[error_level_number] 
                 : d_depth_grad_tol[size-1] );
         size = d_depth_grad_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_depth_grad_time_min[error_level_number] 
                 : d_depth_grad_time_min[size-1] );
         size = d_depth_grad_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_depth_grad_time_max[error_level_number] 
                 : d_depth_grad_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max >= regrid_time);
      }

         //tag based on BATHY_GRADIENT
         if (ref == "BATHY_GRADIENT") {
         var = patch.getPatchData(d_bathy, getDataContext());
         size = d_bathy_grad_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_bathy_grad_tol[error_level_number] 
                 : d_bathy_grad_tol[size-1] );
         size = d_bathy_grad_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_bathy_grad_time_min[error_level_number] 
                 : d_bathy_grad_time_min[size-1] );
         size = d_bathy_grad_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_bathy_grad_time_max[error_level_number] 
                 : d_bathy_grad_time_max[size-1] );
		
         time_allowed = (time_min <= regrid_time) && (time_max >= regrid_time);
      }
			 //tag based on BEDLEVEL_GRADIENT
		    if (ref == "BEDLEVEL_GRADIENT") {
		    var = patch.getPatchData(d_bedlevel, getDataContext());
		    size = d_bedlevel_grad_tol.getSize();
		    tol = ( ( error_level_number < size) 
		            ? d_bedlevel_grad_tol[error_level_number] 
		            : d_bedlevel_grad_tol[size-1] );
		    size = d_bedlevel_grad_time_min.getSize();
		    double time_min = ( ( error_level_number < size) 
		            ? d_bedlevel_grad_time_min[error_level_number] 
		            : d_bedlevel_grad_time_min[size-1] );
		    size = d_bedlevel_grad_time_max.getSize();
		    double time_max = ( ( error_level_number < size) 
		            ? d_bedlevel_grad_time_max[error_level_number] 
		            : d_bedlevel_grad_time_max[size-1] );
	
		    time_allowed = (time_min <= regrid_time) && (time_max >= regrid_time);
		   
		
		 }

	     if (time_allowed) {

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!var.isNull());
#endif

         hier::IntVector vghost   =  var->getGhostCellWidth();
         hier::IntVector tagghost = tags->getGhostCellWidth();
         hier::IntVector ttghost = temp_tags->getGhostCellWidth();

	     //GRADIENT CALC 
         if (ref == "BATHY_GRADIENT" || ref== "DEPTH_GRADIENT" || ref== "BEDLEVEL_GRADIENT") {
            detectgrad_(regrid_time,ifirst(0),ilast(0),   //FORTRAN
                        ifirst(1),ilast(1),
                        vghost(0),vghost(1),
                        tagghost(0),tagghost(1),
                        ttghost(0),ttghost(1),
                        dx,xlo,xhi,
                        tol,
                        TRUE, FALSE,
                        var->getPointer(),
                        tags->getPointer(),temp_tags->getPointer());
         }

      } // if time allowed

   }  // loop over criteria

   /*
    * Adjust temp_tags from those tags set in Richardson extrapolation.
    */
   if (uses_richardson_extrapolation_too) {
        for (pdat::CellIterator ic(pbox); ic; ic++) {              
           if ((*tags)(ic(),0) == RICHARDSON_ALREADY_TAGGED || 
               (*tags)(ic(),0) == RICHARDSON_NEWLY_TAGGED) {
              (*temp_tags)(ic(),0) = TRUE; 
           }
        }
     }

   /*
    * Update tags
    */
   for (pdat::CellIterator ic(pbox); ic; ic++) {
      (*tags)(ic(),0) = (*temp_tags)(ic(),0);
   }

   //stop the timer
   t_taggradient->stop();
}

//======================================================================
//    tagRichardsonExtrapolationCells
//
//    Tag cells for refinement using Richardson Extrapolation with 
//       criteria defined in the input
//   
//    Operates: patch by patch
//
//    this is a concrete implementation of a function from
//    the algs::HyperbolicPatchStrategy<PDIM> abstract base class
//======================================================================

void swe::tagRichardsonExtrapolationCells(
   hier::Patch& patch, 
   const int error_level_number,
   const tbox::Pointer<hier::VariableContext> coarsened_fine, 
   const tbox::Pointer<hier::VariableContext> advanced_coarse,
   const double regrid_time,
   const double deltat,
   const int error_coarsen_ratio,
   const bool initial_error,
   const int tag_index,
   const bool uses_gradient_detector_too)
{
   (void) initial_error;

   const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom = 
      patch.getPatchGeometry();
   const double* xdomainlo = d_grid_geometry->getXLower();
   const double* xdomainhi = d_grid_geometry->getXUpper();

   hier::Box pbox = patch.getBox();

   tbox::Pointer< pdat::CellData<int> > tags     = patch.getPatchData(tag_index);

   // current options for richardson-based tagging:  DEPTH_RICHARDSON
   for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ncrit++) {

      string ref = d_refinement_criteria[ncrit];
      tbox::Pointer< pdat::CellData<double > > coarsened_fine_var; 
      tbox::Pointer< pdat::CellData<double > > advanced_coarse_var;          
      int size = 0;
      double tol = 0.;
      bool time_allowed = false;

      if (ref == "DEPTH_RICHARDSON") {
	     //point to levels of depth
         coarsened_fine_var = patch.getPatchData(d_depth, coarsened_fine);
         advanced_coarse_var = patch.getPatchData(d_depth, advanced_coarse);
         size = d_depth_rich_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_depth_rich_tol[error_level_number] 
                 : d_depth_rich_tol[size-1] );
         size = d_depth_rich_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_depth_rich_time_min[error_level_number] 
                 : d_depth_rich_time_min[size-1] );
         size = d_depth_rich_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_depth_rich_time_max[error_level_number] 
                 : d_depth_rich_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }


      if (time_allowed) {

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!coarsened_fine_var.isNull());
         TBOX_ASSERT(!advanced_coarse_var.isNull());
#endif

         if (ref == "DEPTH_RICHARDSON") {

             //
             // We tag wherever the global error > specified tolerance.
             // The estimated global error is the
             // local truncation error * the approximate number of steps
             // used in the simulation.  Approximate the number of steps as:
             //
             //       steps = L / (s*deltat)
             // where
             //       L = length of problem domain
             //       s = wave speed
             //       delta t = timestep on current level
             //
             // Compute max wave speed from delta t.  This presumes that
             // deltat was computed as deltat = dx/s_max.  We have deltat
             // and dx, so back out s_max from this.
             //

             const double* dx  = patch_geom->getDx();
             double max_dx = 0.;
             double max_length = 0.;
             for (int idir = 0; idir < PDIM; idir++) {
                max_dx = tbox::MathUtilities<double>::Max(max_dx, dx[idir]);
                double length = xdomainhi[idir] - xdomainlo[idir];
                max_length = 
                   tbox::MathUtilities<double>::Max(max_length, length);
             }
             double max_wave_speed = max_dx / deltat;
             double steps = max_length / (max_wave_speed * deltat);


              //
              // Tag cells where |w_c - w_f| * (r^n -1) * steps
              //
              // where
              //       w_c = soln on coarse level (pressure_crse)
              //       w_f = soln on fine level (pressure_fine)
              //       r   = error coarsen ratio
              //       n   = spatial order of scheme (1st or 2nd depending 
              //             on whether Godunov order is 1st or 2nd/4th)
              //
             int order = 1;
             double r = error_coarsen_ratio;
             double rnminus1 = pow(r,order) - 1;

             double diff = 0.;
             double error = 0.;

  
             for (pdat::CellIterator ic(pbox); ic; ic++) {

                //
                // Compute error norm
                //
                diff = (*advanced_coarse_var)(ic(),0) - 
                       (*coarsened_fine_var)(ic(),0);
                error = 
                   tbox::MathUtilities<double>::Abs(diff) * rnminus1 * steps;

                 //
                 // Tag cell if error > prescribed threshold. Since we are
                 // operating on the actual tag values (not temporary ones)
                 // distinguish here tags that were previously set before
                 // coming into this routine and those that are set here.
                 //     RICHARDSON_ALREADY_TAGGED - tagged before coming
                 //                                 into this method.
                 //     RICHARDSON_NEWLY_TAGGED - newly tagged in this method
                 //
                if (error > tol) {
                   if ((*tags)(ic(),0)) {
                      (*tags)(ic(),0) = RICHARDSON_ALREADY_TAGGED;
                   } else {
                      (*tags)(ic(),0) = RICHARDSON_NEWLY_TAGGED;
                   }
                }
            }

         }

      } // time_allowed

   } // loop over refinement criteria

    //
    // If we are NOT performing gradient detector (i.e. only
    // doing Richardson extrapolation) set tags marked in this method
    // to TRUE and all others false.  Otherwise, leave tags set to the 
    // RICHARDSON_ALREADY_TAGGED and RICHARDSON_NEWLY_TAGGED as we may 
    // use this information in the gradient detector.
    //
   if (!uses_gradient_detector_too) {
      for (pdat::CellIterator ic(pbox); ic; ic++) {
         if ( (*tags)(ic(),0) == RICHARDSON_ALREADY_TAGGED ||
             (*tags)(ic(),0) == RICHARDSON_NEWLY_TAGGED ) {
            (*tags)(ic(),0) = TRUE;
         } else {
            (*tags)(ic(),0) = FALSE;
         }
      }
   }

}


//======================================================================
//        registerVisItDataWriter 
//
// register visit data writer to dump output for VisIt
//
//======================================================================

#ifdef HAVE_HDF5
void swe::registerVisItDataWriter(
   tbox::Pointer<appu::VisItDataWriter > viz_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(viz_writer.isNull()));
#endif
   d_visit_writer = viz_writer;
}
#endif

//======================================================================
//        packDerivedDataIntoDoubleBuffer 
//
// add derived quantities (total energy, v*D) to buffer for Vis 
//
//======================================================================

bool swe::packDerivedDataIntoDoubleBuffer(
   double* dbuffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const string& variable_name,
   int zeta_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((region * patch.getBox()) == region);
#endif

   bool data_on_patch = FALSE;
   
   tbox::Pointer< pdat::CellData<double> > depth  =
      patch.getPatchData(d_depth, d_plot_context);
   tbox::Pointer< pdat::CellData<double> > bedlevel  =
      patch.getPatchData(d_bedlevel, d_plot_context);
   tbox::Pointer< pdat::CellData<double> > bathy  =
      patch.getPatchData(d_bathy, d_plot_context);
   tbox::Pointer< pdat::CellData<double> > veldepth =
      patch.getPatchData(d_veldepth, d_plot_context);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!depth.isNull());
   TBOX_ASSERT(!bedlevel.isNull());
   TBOX_ASSERT(!veldepth.isNull());
   TBOX_ASSERT(depth->getGhostBox() == patch.getBox());
   TBOX_ASSERT(veldepth->getGhostBox() == patch.getBox());
#endif

   const hier::Box& data_box = depth->getGhostBox();
   const int box_w0 = region.numberCells(0);
   const int dat_w0 = data_box.numberCells(0);
   const int box_w1 = region.numberCells(1);

   if (variable_name == "Total Energy") {
      const double *const dval  = depth->getPointer(); 
      const double *const udval = veldepth->getPointer(0); 
      const double *const vdval = veldepth->getPointer(1);

      int buf_b1 = 0;
      int dat_b2 = data_box.offset(region.lower());
      int dat_b1 = dat_b2;
      for (int i1 = 0; i1 < box_w1; i1++) {
         for (int i0 = 0; i0 < box_w0; i0++) {
            int dat_indx = dat_b1+i0;
            double v2norm = pow(udval[dat_indx]/dval[dat_indx], 2.0)
                          + pow(vdval[dat_indx]/dval[dat_indx], 2.0)
            ;
            dbuffer[buf_b1+i0] = 
               dval[dat_indx] * (0.5 * v2norm + d_gravity*(dval[dat_indx]));
         } 
         dat_b1 += dat_w0;
         buf_b1 += box_w0;
      }

      data_on_patch = TRUE;

   } else if (variable_name == "velocity") {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(zeta_id < PDIM);
#endif

      const double *const dval = depth->getPointer();
      const double *const veld = veldepth->getPointer(zeta_id);
      int buf_b1 = 0;
      int dat_b2 = data_box.offset(region.lower());

      int dat_b1 = dat_b2;
      for (int i1 = 0; i1 < box_w1; i1++) {
         for (int i0 = 0; i0 < box_w0; i0++) {
            int dat_indx = dat_b1+i0;
            dbuffer[buf_b1+i0] = veld[dat_indx]/(dval[dat_indx]);
         }
         dat_b1 += dat_w0;
         buf_b1 += box_w0;
      }

      data_on_patch = TRUE;

  } else if (variable_name == "bedlevel") {
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(zeta_id < PDIM);
#endif

    const double *const bval = bedlevel->getPointer();
    int buf_b1 = 0;
    int dat_b2 = data_box.offset(region.lower());

    int dat_b1 = dat_b2;
    for (int i1 = 0; i1 < box_w1; i1++) {
       for (int i0 = 0; i0 < box_w0; i0++) {
          int dat_indx = dat_b1+i0;
          dbuffer[buf_b1+i0] = bval[dat_indx];
       }
       dat_b1 += dat_w0;
       buf_b1 += box_w0;
    }

    data_on_patch = TRUE;

  } else if (variable_name == "u") {
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(zeta_id < PDIM);
#endif

	  const double *const dval  = depth->getPointer();
	  const double *const udval = veldepth->getPointer(0);
	  int buf_b1 = 0;
	  int dat_b2 = data_box.offset(region.lower());

	  int dat_b1 = dat_b2;
	  for (int i1 = 0; i1 < box_w1; i1++) {
	     for (int i0 = 0; i0 < box_w0; i0++) {
	       int dat_indx = dat_b1+i0;
	        dbuffer[buf_b1+i0] = udval[dat_indx]/(dval[dat_indx]);
	     }
         dat_b1 += dat_w0;
	     buf_b1 += box_w0;
	  }

	  data_on_patch = TRUE;
	
	  } else if (variable_name == "v") {
	#ifdef DEBUG_CHECK_ASSERTIONS
	   TBOX_ASSERT(zeta_id < PDIM);
	#endif

		  const double *const dval  = depth->getPointer();
		  const double *const udval = veldepth->getPointer(1);
		  int buf_b1 = 0;
		  int dat_b2 = data_box.offset(region.lower());

		  int dat_b1 = dat_b2;
		  for (int i1 = 0; i1 < box_w1; i1++) {
		     for (int i0 = 0; i0 < box_w0; i0++) {
		       int dat_indx = dat_b1+i0;
		        dbuffer[buf_b1+i0] = udval[dat_indx]/(dval[dat_indx]);
		     }
	         dat_b1 += dat_w0;
		     buf_b1 += box_w0;
		  }

		  data_on_patch = TRUE;
		
		 } else if (variable_name == "zeta") {
		#ifdef DEBUG_CHECK_ASSERTIONS
		   TBOX_ASSERT(zeta_id < PDIM);
		#endif

			  const double *const dval = depth->getPointer();
			  const double *const bval = bathy->getPointer();
			  int buf_b1 = 0;
			  int dat_b2 = data_box.offset(region.lower());

			  int dat_b1 = dat_b2;
			  for (int i1 = 0; i1 < box_w1; i1++) {
			     for (int i0 = 0; i0 < box_w0; i0++) {
			       int dat_indx = dat_b1+i0;
			        dbuffer[buf_b1+i0] = bval[dat_indx]+dval[dat_indx];
			     }
		         dat_b1 += dat_w0;
			     buf_b1 += box_w0;
			  }

			  data_on_patch = TRUE;
			
		} else if (variable_name == "wetdry") {

			 const double *const dval = depth->getPointer();
			 int buf_b1 = 0;
			 int dat_b2 = data_box.offset(region.lower());

			 int dat_b1 = dat_b2;
			 for (int i1 = 0; i1 < box_w1; i1++) {
				for (int i0 = 0; i0 < box_w0; i0++) {
				   int dat_indx = dat_b1+i0;
			   	   if(dval[dat_indx] < 1e-6)
			         {
				     dbuffer[buf_b1+i0] = 0.0;
			         }
				   else
				     {
				     dbuffer[buf_b1+i0] = 1.0;
			         }
				 }
			     dat_b1 += dat_w0;
				 buf_b1 += box_w0;
			  }

			  data_on_patch = TRUE;
			


   } else {
      TBOX_ERROR("swe::packDerivedDataIntoDoubleBuffer()"
             << "\n    unknown variable_name " << variable_name << "\n");
   }

   return(data_on_patch);
   
}

//======================================================================
//        writeData1dPencil     
//
// write 1d data on a line through patch along a 'pencil box' to file 
//
//======================================================================

void swe::writeData1dPencil(const tbox::Pointer<hier::Patch > patch,
                              const hier::Box& pencil_box,
                              const int idir,
                              ostream& file)
{

   const hier::Box& patchbox = patch->getBox();
   const hier::Box box = pencil_box * patchbox;

   if (!box.empty()) {

      tbox::Pointer< pdat::CellData<double> > depth  =
         patch->getPatchData(d_depth, getDataContext());
      tbox::Pointer< pdat::CellData<double> > bathy  =
         patch->getPatchData(d_bathy, getDataContext());
      tbox::Pointer< pdat::CellData<double> > veldepth =
         patch->getPatchData(d_veldepth, getDataContext());
      tbox::Pointer< pdat::CellData<double> > bedlevel =
         patch->getPatchData(d_bedlevel, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!depth.isNull());
      TBOX_ASSERT(!veldepth.isNull());
      TBOX_ASSERT(!bedlevel.isNull());
#endif

      const tbox::Pointer<geom::CartesianPatchGeometry > pgeom = patch->getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();

      const double cell_center = xlo[idir] 
                                 + (double(box.lower(idir)
                                     - patchbox.lower(idir)) 
                                    + 0.5) * dx[idir];

      int ccount = 0;
      for (pdat::CellIterator ic(box); ic; ic++) {
               file << cell_center + ccount*dx[idir] << " ";
               ccount++;
               double dval  = (*depth)(ic(),0);
			   double bval  = (*bathy)(ic(),0);
               double udval = (*veldepth)(ic(),idir);
               double vel  = udval/dval;
			   double wetdry = 0; // wet
				double blval = (*bedlevel)(ic(),0);
		       if(dval < 1e-30){
				wetdry = 1; 
				vel = 0;
			}
               //double etot = d_rho*(0.5*vel*vel + d_gravity*dval);
           
               /*
                * Write out conserved quantities.
                */ 
               file << dval << " ";
           	   file << vel << " ";
			   file << dval+bval << " ";
			   file << bval << " ";
			   file << wetdry << " ";
				file << blval << " ";
               file << endl;
            }

   }

}

//======================================================================
//        printClassData        
//
//    dump general info (parameters and such) to output 
//
//======================================================================

void swe::printClassData(ostream &os) const 
{
   int j;

   os << "\nswe::printClassData..." << endl;
   os << "swe: this = " << (swe*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometry = "
      << (geom::CartesianGridGeometry *)d_grid_geometry << endl;

 
   os << "Ghost Sizes..." << endl;
   os << "   d_nghosts = " << d_nghosts << endl;
   os << "   d_fluxghosts = " << d_fluxghosts << endl;

   os << "Problem description and initial data..." << endl;
   os << "   d_data_problem = " << d_data_problem << endl;
   os << "   d_data_problem_int = " << d_data_problem_int << endl;
   os << "   d_C_manning  =  " << d_C_manning << endl;
   os << "   d_mindepth   =  " << d_mindepth  << endl;
   os << "   d_frictype   =  " << d_frictype  << endl;
   os << "   d_fluxorder  =  " << d_fluxorder << endl;
   os << "   d_transverse =  " << d_transverse << endl;
   os << "   d_sedmodel   =  " << d_sedmodel << endl; 
   os << "   d_sedinit    =  " << d_sedinit << endl;
   os << "   d_taucrit    =  " << d_taucrit << endl;
   os << "   d_morphfactor  =  " << d_morphfactor << endl;



   os << "   Boundary condition data " << endl;
   os << " Dirichlet BC is type " << DIRICHLET_BC << endl ;
   os << " Reflect BC is type " << REFLECT_BC << endl ;
   os << " Neumann BC is type " << NEUMANN_BC << endl ;
   os << " FLOW BC is type " << FLOW_BC << endl ;

   for (j = 0; j < d_master_bdry_edge_conds.getSize(); j++) {
      os << "\n       d_master_bdry_edge_conds[" << j << "] = "
         << d_master_bdry_edge_conds[j] << endl;
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      os << "       d_vector_bdry_edge_conds[" << j << "] = "
         << d_vector_bdry_edge_conds[j] << endl;
      if (d_master_bdry_edge_conds[j] == DIRICHLET_BC) {
         os << "         d_bdry_edge_depth[" << j << "] = "
            << d_bdry_edge_depth[j] << endl;
         os << "         d_bdry_edge_veldepth[" << j << "] = "
            << d_bdry_edge_veldepth[j*PDIM+0] << " , "
            << d_bdry_edge_veldepth[j*PDIM+1] << endl;
      }
   }
   os << endl;
   for (j = 0; j < d_master_bdry_node_conds.getSize(); j++) {
      os << "\n       d_master_bdry_node_conds[" << j << "] = "
         << d_master_bdry_node_conds[j] << endl;
      os << "       d_scalar_bdry_node_conds[" << j << "] = "
         << d_scalar_bdry_node_conds[j] << endl;
      os << "       d_vector_bdry_node_conds[" << j << "] = "
         << d_vector_bdry_node_conds[j] << endl;
      os << "       d_node_bdry_edge[" << j << "] = "
         << d_node_bdry_edge[j] << endl;
   }

   os << "   Refinement criteria parameters " << endl;
 
   for (j = 0; j < d_refinement_criteria.getSize(); j++) {
      os << "       d_refinement_criteria[" << j << "] = "
         << d_refinement_criteria[j] << endl;
   }
   for (j = 0; j < d_depth_grad_tol.getSize(); j++) {
      os << "       d_depth_grad_tol[" << j << "] = "
         << d_depth_grad_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_depth_grad_time_max.getSize(); j++) {
      os << "       d_depth_grad_time_max[" << j << "] = "
         << d_depth_grad_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_depth_grad_time_min.getSize(); j++) {
      os << "       d_depth_grad_time_min[" << j << "] = "
         << d_depth_grad_time_min[j] << endl;
   }
   for (j = 0; j < d_bedlevel_grad_tol.getSize(); j++) {
      os << "       d_bedlevel_grad_tol[" << j << "] = "
         << d_bedlevel_grad_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_bedlevel_grad_time_max.getSize(); j++) {
      os << "       d_bedlevel_grad_time_max[" << j << "] = "
         << d_bedlevel_grad_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_bedlevel_grad_time_min.getSize(); j++) {
      os << "       d_bedlevel_grad_time_min[" << j << "] = "
         << d_bedlevel_grad_time_min[j] << endl;
   }
  
   os << endl;
   for (j = 0; j < d_depth_rich_tol.getSize(); j++) {
      os << "       d_depth_rich_tol[" << j << "] = "
         << d_depth_rich_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_depth_rich_time_max.getSize(); j++) {
      os << "       d_depth_rich_time_max[" << j << "] = "
         << d_depth_rich_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_depth_rich_time_min.getSize(); j++) {
      os << "       d_depth_rich_time_min[" << j << "] = "
         << d_depth_rich_time_min[j] << endl;
   }
   os << endl;

}

//======================================================================
//        getFromInput          
//
//    get runtime parameters from input 
//
//    I believe this overrides the vals from restart so be careful 
//
//    this is a concrete implementation of a function from
//    the ? class 
//
//======================================================================
void swe::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   //from Euler, if restarting, nonuniform workload only allowed if
   //used prior
   if (!is_from_restart) { 
      d_use_nonuniform_workload = 
         db->getBoolWithDefault("use_nonuniform_workload",
                                d_use_nonuniform_workload);
   } else {
      if (d_use_nonuniform_workload) {
         d_use_nonuniform_workload = 
            db->getBool("use_nonuniform_workload");
      }
   }

 
   //refinement data options
   if (db->keyExists("Refinement_data")) {
      tbox::Pointer<tbox::Database> refine_db = db->getDatabase("Refinement_data");
      tbox::Array<string> refinement_keys = refine_db->getAllKeys();
      int num_keys = refinement_keys.getSize();

      if (refine_db->keyExists("refine_criteria")) {
         d_refinement_criteria = 
            refine_db->getStringArray("refine_criteria");
      } else {
         TBOX_WARNING(d_object_name << ": "
                   << "No key `refine_criteria' found in data for"
                   << " RefinementData. No refinement will occur." << endl);
      }
       
      tbox::Array<string> ref_keys_defined(num_keys);
      int def_key_cnt = 0;
      tbox::Pointer<tbox::Database> error_db;
      for (int i = 0; i < refinement_keys.getSize(); i++) {
           
         string error_key = refinement_keys[i];
         error_db.setNull();

         if ( !(error_key == "refine_criteria") ) {

            if ( !(error_key == "DEPTH_GRADIENT" || error_key=="BATHY_GRADIENT" || error_key=="BEDLEVEL_GRADIENT") ){
               TBOX_ERROR(d_object_name << ": "
                         << "Unknown refinement criteria: " << error_key
                         << "\nin input." << endl);
            } else {
               error_db = refine_db->getDatabase(error_key);
               ref_keys_defined[def_key_cnt] = error_key;
               def_key_cnt++;
            }
               
              
            if (!error_db.isNull() && error_key == "DEPTH_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_depth_grad_tol = 
                  error_db->getDoubleArray("grad_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `grad_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_depth_grad_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_depth_grad_time_max.resizeArray(1);
                  d_depth_grad_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_depth_grad_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_depth_grad_time_min.resizeArray(1);
                  d_depth_grad_time_min[0] = 0.;
               } 

            }
            if (!error_db.isNull() && error_key == "BATHY_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_bathy_grad_tol = 
                  error_db->getDoubleArray("grad_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `grad_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_bathy_grad_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_bathy_grad_time_max.resizeArray(1);
                  d_bathy_grad_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_bathy_grad_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_bathy_grad_time_min.resizeArray(1);
                  d_bathy_grad_time_min[0] = 0.;
               } 

            }
				 if (!error_db.isNull() && error_key == "BEDLEVEL_GRADIENT") {

		               if (error_db->keyExists("grad_tol")) {
		                  d_bedlevel_grad_tol = 
		                  error_db->getDoubleArray("grad_tol");
		               } else {
		                  TBOX_ERROR(d_object_name << ": "
		                           << "No key `grad_tol' found in data for "
		                           << error_key << endl);
		               }

		               if (error_db->keyExists("time_max")){
		                  d_bedlevel_grad_time_max = 
		                  error_db->getDoubleArray("time_max");
		               } else {
		                  d_bedlevel_grad_time_max.resizeArray(1);
		                  d_bedlevel_grad_time_max[0] = 
		                     tbox::MathUtilities<double>::getMax();
		               }

		               if (error_db->keyExists("time_min")){
		                  d_bedlevel_grad_time_min = 
		                  error_db->getDoubleArray("time_min");
		               } else {
		                  d_bedlevel_grad_time_min.resizeArray(1);
		                  d_bedlevel_grad_time_min[0] = 0.;
		               } 

		            }
              

            if (!error_db.isNull() && error_key == "DEPTH_RICHARDSON") {

               if (error_db->keyExists("rich_tol")) {
                  d_depth_rich_tol = 
                  error_db->getDoubleArray("rich_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `rich_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_depth_rich_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_depth_rich_time_max.resizeArray(1);
                  d_depth_rich_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_depth_rich_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_depth_rich_time_min.resizeArray(1);
                  d_depth_rich_time_min[0] = 0.;
               } 

            }


         } 

      } //loop over refinement_keys.getSize

      for (int k0 = 0; k0 < d_refinement_criteria.getSize(); k0++) {
         string use_key = d_refinement_criteria[k0];
         bool key_found = false;
         for (int k1 = 0; k1 < def_key_cnt; k1++) {
             string def_key = ref_keys_defined[k1];
             if (def_key == use_key) key_found = true;
         }

         if (!key_found) {
             TBOX_ERROR(d_object_name << ": "
                       << "No input found for specified refine criteria: "
                       << d_refinement_criteria[k0] << endl);
         }
      }

   } // if "Refinement_data" db entry exists

   //specify the problem to solve, a potential intial data database
   //do not pull from input if we are restarting (can't change problem on the fly)
   if (!is_from_restart) { 

      if (db->keyExists("data_problem")) {
         d_data_problem = db->getString("data_problem");
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`data_problem' value not found in input." << endl);
      }

      // Read the Manning coefficient.  
      if (db->keyExists("C_manning")) {
         d_C_manning = db->getDouble("C_manning"); 
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`C_manning' value not found in input." << endl);
      }

      // Read the minimum depth 
      if (db->keyExists("mindepth")) {
         d_mindepth = db->getDouble("mindepth");
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`mindepth' value not found in input." << endl);
      }

      // Read the Friction Type 
      if (db->keyExists("frictype")) {
         d_frictype  = db->getInteger("frictype");
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`frictype' value not found in input." << endl);
      }

      // Read the flux order    
      if (db->keyExists("fluxorder")) {
         d_fluxorder  = db->getInteger("fluxorder");
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`fluxorder' value not found in input." << endl);
      }

      // Read the transverse correction control var
      if (db->keyExists("transverse")) {
         d_transverse = db->getInteger("transverse");
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`transverse' value not found in input." << endl);
      }

     // Read the sediment model activation var
      if (db->keyExists("sedmodel")) {
         d_sedmodel = db->getInteger("sedmodel");
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`sedmodel' value not found in input." << endl);
      }

		// Read the sediment model activation var
		if (db->keyExists("taucrit")) {
		   d_taucrit = db->getDouble("taucrit");
		} else {
		   TBOX_ERROR(d_object_name << ": "
		      << "`taucrit' value not found in input." << endl);
		}

		// Read the sediment model initial Time var
		if (db->keyExists("sedinit")) {
	      d_sedinit = db->getDouble("sedinit");
	   } else {
	      TBOX_ERROR(d_object_name << ": "
	         << "`sedinit' value not found in input." << endl);
	   }
		
     // Read the sediment model activation var
      if (db->keyExists("morphfactor")) {
         d_morphfactor = db->getDouble("morphfactor");
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`morphfactor' value not found in input." << endl);
      }









      //tbox::Pointer<tbox::Database> init_data_db;
      //if (db->keyExists("Initial_data")) {
      //   init_data_db = db->getDatabase("Initial_data");
      //} else {
      //   TBOX_ERROR(d_object_name << ": "
      //            << "No `Initial_data' database found in input." << endl);
      //}

      bool found_problem_data = false;

      //here is where we would pull data params for a specific problem
      //if (d_data_problem == "SPHERE") {
      //
      //   if (init_data_db->keyExists("radius")) {
      //      d_radius = init_data_db->getDouble("radius");
      //   } else {
      //      TBOX_ERROR(d_object_name << ": "
      //         << "`radius' input required for SPHERE problem." << endl);
      //   }
   
      //right now we don't have specific problem data
      found_problem_data = true;

      if (!found_problem_data) {
         TBOX_ERROR(d_object_name << ": "
            << "`Initial_data' database found in input." 
            << " But bad data supplied." << endl);
      } 
 
   } // if !is_from_restart read in problem data

   //read the boundary data from input
   const hier::IntVector& one_vec = hier::IntVector::getOne(d_dim);
   hier::IntVector periodic = d_grid_geometry->getPeriodicShift(one_vec);
   int num_per_dirs = 0;
   for (int id = 0; id < PDIM; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < PDIM) {

      if (db->keyExists("Boundary_data")) {

         tbox::Pointer<tbox::Database> bdry_db = db->getDatabase("Boundary_data");

			if(d_dim == tbox::Dimension(2)){
         appu::CartesianBoundaryUtilities2::readBoundaryInput(this,
                                                        bdry_db,
                                                        d_master_bdry_edge_conds,
                                                        d_master_bdry_node_conds,
                                                        periodic);
			}

      } else {
         TBOX_ERROR(d_object_name << ": "
            << "Key data `Boundary_data' not found in input. " << endl);
      }

   }

}

//======================================================================
//        putToDataBase         
//
//    put runtime parameters to restart database
//
//    this is a concrete implementation of a function from
//    the ? class 
//======================================================================


void swe::putToDatabase(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("SWE_VERSION", SWE_VERSION);

   db->putIntegerArray("d_nghosts", &d_nghosts[0], PDIM);
   db->putIntegerArray("d_fluxghosts", &d_fluxghosts[0], PDIM);
   db->putString("d_data_problem", d_data_problem);
   db->putDouble("d_C_manning", d_C_manning);
   db->putDouble("d_mindepth",  d_mindepth);  
   db->putInteger("d_fluxorder",d_fluxorder);
   db->putInteger("d_transverse",d_transverse);
   db->putInteger("d_frictype",d_frictype);
   db->putInteger("d_sedmodel",d_sedmodel);
   db->putDouble("d_sedinit",d_sedinit);
   db->putDouble("d_taucrit",d_taucrit);
   db->putDouble("d_morphfactor",d_morphfactor);

   db->putIntegerArray("d_master_bdry_edge_conds", d_master_bdry_edge_conds);
   db->putIntegerArray("d_master_bdry_node_conds", d_master_bdry_node_conds);

   db->putDoubleArray("d_bdry_edge_depth", d_bdry_edge_depth);
   db->putDoubleArray("d_bdry_edge_bedlevel", d_bdry_edge_bedlevel);
   db->putDoubleArray("d_bdry_edge_bathy", d_bdry_edge_bathy);
   db->putDoubleArray("d_bdry_edge_veldepth", d_bdry_edge_veldepth);

   if (d_refinement_criteria.getSize() > 0) {
      db->putStringArray("d_refinement_criteria", d_refinement_criteria);
   }
   for (int i = 0; i < d_refinement_criteria.getSize(); i++) {


      if (d_refinement_criteria[i] == "DEPTH_GRADIENT") {
		 db->putDoubleArray("d_depth_grad_tol", 
                            d_depth_grad_tol);
         db->putDoubleArray("d_depth_grad_time_max", 
                            d_depth_grad_time_max);
         db->putDoubleArray("d_depth_grad_time_max", 
                            d_depth_grad_time_max);

      } else if  (d_refinement_criteria[i] == "DEPTH_RICHARDSON") {

         db->putDoubleArray("d_depth_rich_tol", 
                            d_depth_rich_tol);
         db->putDoubleArray("d_depth_rich_time_max", 
                            d_depth_rich_time_max);
         db->putDoubleArray("d_depth_rich_time_min", 
                            d_depth_rich_time_min);
 
      }
   }

}

//======================================================================
//        getFromRestart        
//
//    get runtime parameters to restart database
//
//    this is a concrete implementation of a function from
//    the ? class 
//======================================================================
void swe::getFromRestart()
{

   tbox::Pointer<tbox::Database> root_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
                 << d_object_name << " not found in restart file." << endl);
   }

   int ver = db->getInteger("SWE_VERSION");
   if (ver != SWE_VERSION) {
      TBOX_ERROR(d_object_name << ": "
          << "Restart file version different than class version." << endl);
   }

 
   int* tmp_nghosts = &d_nghosts[0];
   db->getIntegerArray("d_nghosts", tmp_nghosts, PDIM);
   for (int i = 0; i < PDIM; i++) {
      if (d_nghosts(i) != CELLG) {
         TBOX_ERROR(d_object_name << ": "
            << "Key data `d_nghosts' in restart file != CELLG." << endl);
      }
   }
   int* tmp_fluxghosts = &d_fluxghosts[0];
   db->getIntegerArray("d_fluxghosts", tmp_fluxghosts, PDIM);
   for (int i = 0; i < PDIM; i++) {
      if (d_fluxghosts(i) != FLUXG) {
         TBOX_ERROR(d_object_name << ": "
            << "Key data `d_fluxghosts' in restart file != FLUXG." << endl);
      }
   }

   d_data_problem = db->getString("d_data_problem");
   d_C_manning    = db->getDouble("d_C_manning");
   d_mindepth     = db->getDouble("d_mindepth");
   d_fluxorder    = db->getInteger("d_fluxorder");
   d_frictype     = db->getInteger("d_frictype");
   d_frictype     = db->getInteger("d_transverse");
   d_sedmodel     = db->getInteger("d_sedmodel");
   d_sedinit      = db->getDouble("d_sedinit");
   d_taucrit      = db->getDouble("d_taucrit");
   d_morphfactor  = db->getDouble("d_morphfactor");

   d_master_bdry_edge_conds = db->getIntegerArray("d_master_bdry_edge_conds");
   d_master_bdry_node_conds = db->getIntegerArray("d_master_bdry_node_conds");

   d_bdry_edge_depth = db->getDoubleArray("d_bdry_edge_depth");
   d_bdry_edge_bedlevel = db->getDoubleArray("d_bdry_edge_bedlevel");
   d_bdry_edge_bathy = db->getDoubleArray("d_bdry_edge_bathy");
   d_bdry_edge_veldepth = db->getDoubleArray("d_bdry_edge_veldepth");

   if (db->keyExists("d_refinement_criteria")) {
      d_refinement_criteria = db->getStringArray("d_refinement_criteria");
   }

   for (int i = 0; i < d_refinement_criteria.getSize(); i++) {

      if (d_refinement_criteria[i] == "DEPTH_GRADIENT") {

         d_depth_grad_tol = 
            db->getDoubleArray("d_depth_grad_tol");
         d_depth_grad_time_max = 
            db->getDoubleArray("d_depth_grad_time_max");
         d_depth_grad_time_min = 
            db->getDoubleArray("d_depth_grad_time_min");

      } else if  (d_refinement_criteria[i] == "depth_RICHARDSON") {

         d_depth_rich_tol = 
            db->getDoubleArray("d_depth_rich_tol");
         d_depth_rich_time_max = 
            db->getDoubleArray("d_depth_rich_time_max");
         d_depth_rich_time_min = 
            db->getDoubleArray("d_depth_rich_time_min");

      }

   }

}

//======================================================================
//        readDirichletBoundaryDataEntry       
//
//    read boundary data from input database
//
//======================================================================

void swe::readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                           string& db_name,
                                           int bdry_location_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
   TBOX_ASSERT(!db_name.empty());
#endif

	if (d_dim == tbox::Dimension(2)) {
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_edge_depth,
                      d_bdry_edge_bathy,
                      d_bdry_edge_bedlevel,
					  		 d_bdry_edge_veldepth);
	}
      
}

void swe::readStateDataEntry(tbox::Pointer<tbox::Database> db,
                               const string& db_name,
                               int array_indx,
                               tbox::Array<double>& depth,
                               tbox::Array<double>& bathy,
										 tbox::Array<double>& bedlevel,
                               tbox::Array<double>& veldepth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
   TBOX_ASSERT(!db_name.empty());
   TBOX_ASSERT(array_indx >= 0);
   TBOX_ASSERT(depth.getSize() > array_indx);
   TBOX_ASSERT(bathy.getSize() > array_indx);
   TBOX_ASSERT(bedlevel.getSize() > array_indx);
   TBOX_ASSERT(veldepth.getSize() > array_indx*PDIM);
#endif

   if (db->keyExists("depth")) {
      depth[array_indx] = db->getDouble("depth");
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`depth' entry missing from " << db_name
         << " input database. " << endl);
   }
   if (db->keyExists("bathy")) {
      bathy[array_indx] = db->getDouble("bathy");
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`bathy' entry missing from " << db_name
         << " input database. " << endl);
   }
   if (db->keyExists("bedlevel")) {
      bedlevel[array_indx] = db->getDouble("bedlevel");
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`bedlevel' entry missing from " << db_name
         << " input database. " << endl);
   }
   if (db->keyExists("veldepth")) {
      tbox::Array<double> tmp_vel(0);
      tmp_vel = db->getDoubleArray("veldepth");
      if (tmp_vel.getSize() < PDIM) {
         TBOX_ERROR(d_object_name << ": "
                    << "Insufficient number `veldepth' values"
                    << " given in " << db_name
                    << " input database." << endl);
      }
      for (int iv = 0; iv < PDIM; iv++) {
         veldepth[array_indx*PDIM+iv] = tmp_vel[iv];
      }
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`veldepth' entry missing from " << db_name
         << " input database. " << endl);
   }
  

}

//======================================================================
//        checkBoundaryData     
//
//    check the boundary data when debugging
//
//======================================================================


void swe::checkBoundaryData(int btype,
                              const hier::Patch& patch,
                              const hier::IntVector& ghost_width_to_check,
                              const tbox::Array<int>& scalar_bconds,
                              const tbox::Array<int>& vector_bconds) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(btype == EDGE2D_BDRY_TYPE || 
          btype == NODE2D_BDRY_TYPE);
#endif

   const tbox::Pointer<geom::CartesianPatchGeometry > pgeom = patch.getPatchGeometry();
   const tbox::Array<hier::BoundaryBox> bdry_boxes =
      pgeom->getCodimensionBoundaries(btype);

   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

   for (int i = 0; i < bdry_boxes.getSize(); i++ ) {
      hier::BoundaryBox bbox = bdry_boxes[i];
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(bbox.getBoundaryType() == btype);
#endif
      int bloc = bbox.getLocationIndex();

      int bscalarcase, bveldepthcase, refbdryloc;

		if (d_dim == tbox::Dimension(2)) {
			if (btype == EDGE2D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
				TBOX_ASSERT(scalar_bconds.getSize() == NUM_2D_EDGES);
				TBOX_ASSERT(vector_bconds.getSize() == NUM_2D_EDGES);
#endif
				bscalarcase = scalar_bconds[bloc];
				bveldepthcase = vector_bconds[bloc];
				refbdryloc = bloc;
			} else { // btype == NODE2D_BDRY_TYPE
#ifdef DEBUG_CHECK_ASSERTIONS
				TBOX_ASSERT(scalar_bconds.getSize() == NUM_2D_NODES);
				TBOX_ASSERT(vector_bconds.getSize() == NUM_2D_NODES);
#endif
				bscalarcase = scalar_bconds[bloc];
				bveldepthcase = vector_bconds[bloc];
				refbdryloc = d_node_bdry_edge[bloc];
			}
		}


      int num_bad_values = 0;

      num_bad_values = 
      appu::CartesianBoundaryUtilities2::checkBdryData(
         d_depth->getName(),
         patch,
         vdb->mapVariableAndContextToIndex(d_depth, getDataContext()), 0, 
         ghost_width_to_check,
         bbox,
         bscalarcase,
         d_bdry_edge_depth[refbdryloc]);
#if (TESTING == 1)
      if (num_bad_values > 0) {
         tbox::perr << "\nswe Boundary Test FAILED: \n" 
           << "     " << num_bad_values << " bad depth values found for\n"
           << "     boundary type " << btype << " at location " << bloc << endl;
      }
#endif





      for (int idir = 0; idir < PDIM; idir++) {

         int vbcase = bscalarcase;
			if (d_dim == tbox::Dimension(2)) {
	         if (btype == EDGE2D_BDRY_TYPE) {
            	if ( (idir == 0 && (bloc == XLO || bloc == XHI)) ||
                 	(idir == 1 && (bloc == YLO || bloc == YHI)) ) {
               	vbcase = bveldepthcase;
            	}
				} else if (btype == NODE2D_BDRY_TYPE) {            
					if ( (idir == 0 && bveldepthcase == XREFLECT_BC) ||
					(idir == 1 && bveldepthcase == YREFLECT_BC) ) {
					vbcase = bveldepthcase;
					}
				}
			}


			if (d_dim == tbox::Dimension(2)) {
         	num_bad_values =
         	appu::CartesianBoundaryUtilities2::checkBdryData(
               d_veldepth->getName(),
               patch,
               vdb->mapVariableAndContextToIndex(d_veldepth, getDataContext()), 
               idir,
               ghost_width_to_check,
               bbox,
               vbcase,
               d_bdry_edge_veldepth[refbdryloc*PDIM+idir]);
			}

#if (TESTING == 1)
         if (num_bad_values > 0) {
            tbox::perr << "\nswe Boundary Test FAILED: \n" 
              << "     " << num_bad_values 
              << " bad VELDEPTH values found in direction " << idir << " for\n"
              << "     boundary type " << btype << " at location " << bloc << endl;
         }
#endif
      }

   }

}
