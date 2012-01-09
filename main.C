
// 
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/examples/swe/main.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Main program for SAMRAI swe gas dynamics sample application
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
using namespace std;

#ifdef OPENCL
#include <CL/cl.h>
#endif 

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <sys/stat.h>

// Headers for major algorithm/data structure objects from SAMRAI

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/algs/TimeRefinementIntegrator.h"
#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"

// Headers for basic SAMRAI objects

#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/SiloDatabaseFactory.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"


// Header for application-specific algorithm/data structure object

#include "swe.h"

// Classes for autotesting.

#if (TESTING == 1)
#include "AutoTester.h"
#endif

using namespace SAMRAI;

/*
 ************************************************************************
 *                                                                      *
 * This is the main program for an AMR shallow water eqns application   *
 * built using SAMRAI.   The application program is constructed by      *
 * composing a variety of algorithm objects found in SAMRAI plus some   *
 * others that are specific to this application. The following brief    *
 * discussion summarizes these objects.                                 *
 *                                                                      *
 *    hier::PatchHierarchy - A container for the AMR patch hierarchy and      *
 *       the data on the grid.                                          *
 *                                                                      *
 *    geom::CartesianGridGeometry<NDIM> - Defines and maintains the Cartesian       *
 *       coordinate system on the grid.  The hier::PatchHierarchy<NDIM>             *
 *       maintains a reference to this object.                          *
 *                                                                      *
 * A single overarching algorithm object drives the time integration    *
 * and adaptive gridding processes:                                     *
 *                                                                      *
 *    algs::TimeRefinementIntegrator<NDIM> - Coordinates time integration and       *
 *       adaptive gridding procedures for the various levels            *
 *       in the AMR patch hierarchy.  Local time refinement is          *
 *       employed during hierarchy integration; i.e., finer             *
 *       levels are advanced using smaller time increments than         *
 *       coarser level.  Thus, this object also invokes data            *
 *       synchronization procedures which couple the solution on        *
 *       different patch hierarchy levels.                              *
 *                                                                      *
 * The time refinement integrator is not specific to the numerical      *
 * methods used and the problem being solved.   It maintains references *
 * to two other finer grain algorithmic objects, more specific to       *
 * the problem at hand, with which it is configured when they are       *
 * passed into its constructor.   They are:                             *
 *                                                                      *
 *    algs::HyperbolicLevelIntegrator<NDIM> - Defines data management procedures    *
 *       for level integration, data synchronization between levels,    *
 *       and tagging cells for refinement.  These operations are        *
 *       tailored to explicit time integration algorithms used for      *
 *       hyperbolic systems of conservation laws, such as the shallow   *
 *       water eqns. This integrator manages data for numerical         *
 *       routines that treat individual patches in the AMR patch        *
 *       hierarchy.  In this particular application, it maintains a     *
 *       pointer to the swe object that defines variables and           *
 *       provides numerical routines for the swe model.                 *
 *                                                                      *
 *       swe- Defines variables and numerical routines for the          *
 *          discrete shallow water eqns on each patch in the AMR        *
 *          hierarchy.                                                  *
 *                                                                      *
 *    mesh::GriddingAlgorithm<NDIM> - Drives the AMR patch hierarchy generation     *
 *       and regridding procedures.  This object maintains              *
 *       references to three other algorithmic objects with             *
 *       which it is configured when they are passed into its           *
 *       constructor.   They are:                                       *
 *                                                                      *
 *       mesh::BergerRigoutsos<NDIM> - Clusters cells tagged for refinement on a    *
 *          patch level into a collection of logically-rectangular      *
 *          box domains.                                                *
 *                                                                      *
 *       mesh::LoadBalancer<NDIM> - Processes the boxes generated by the            *
 *          mesh::BergerRigoutsos<NDIM> algorithm into a configuration from         *
 *          which patches are contructed.  The algorithm we use in this *
 *          class assumes a spatially-uniform workload distribution;    *
 *          thus, it attempts to produce a collection of boxes          *
 *          each of which contains the same number of cells.  The       *
 *          load balancer also assigns patches to processors.           *
 *                                                                      *
 *       mesh::StandardTagAndInitialize<NDIM> - Couples the gridding algorithm      *
 *          to the HyperbolicIntegrator. Selects cells for              *
 *          refinement based on either Gradient detection, Richardson   *
 *          extrapolation, or pre-defined Refine box region.  The       *
 *          object maintains a pointer to the algs::HyperbolicLevelIntegrator<NDIM>,*
 *          which is passed into its constructor, for this purpose.     *
 *                                                                      *
 ************************************************************************
 */

/* 
 *******************************************************************
 *                                                                 *
 * For each run, the input filename and restart information        *
 * (if needed) must be given on the command line.                  *
 *                                                                 *
 *      For non-restarted case, command line is:                   *
 *                                                                 *
 *          executable <input file name>                           *
 *                                                                 *
 *      For restarted run, command line is:                        *
 *                                                                 *
 *          executable <input file name> <restart directory> \     *
 *                     <restart number>                            *
 *                                                                 *
 * Accessory routines used within the main program:                *
 *                                                                 *
 *   dumpVizData1dPencil - Writes 1d pencil of swe solution data   *
 *      to plot files so that it may be viewed in MatLab.  This    *
 *      routine assumes a single patch level in 2d and 3d.  In     *
 *      other words, it only plots data on level zero.  It can     *
 *      handle AMR in 1d.                                          *
 *                                                                 *
 *******************************************************************
 */

static void dumpMatlabData1dPencil(const string& dirname, 
                            const string& filename,
                            const int ext,
                            const double plot_time,
                            const tbox::Pointer<hier::PatchHierarchy > hierarchy,
                            const int pencil_direction,
                            const bool default_pencil,
                            const tbox::Array<int>& pencil_index,
                            swe* swe_model);

						
void dumpProbe(const string& dirname, 
						    const string& filename,
						    const int ext,
						    const double plot_time,
						    const tbox::Pointer<hier::PatchHierarchy > hierarchy,
						    const double probe_loc_x,
						    const double probe_loc_y,
						    swe* swe_model);

//======================================================================
//          main program 
//======================================================================

int main( int argc, char *argv[] )
{
	

   //initialize MPI and SAMRAI
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   int num_failures = 0;

   {

   string input_filename;
   string restart_read_dirname;
   int restore_num = 0;

   bool is_from_restart = false;
	 bool uses_visit = true;

   //process commandline
   if ( (argc != 2) && (argc != 4) ) {
        tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
           << "<restart dir> <restore number> [options]\n"
           << "  options:\n"
           << "  none at this time"
           << endl;
        tbox::SAMRAI_MPI::abort();
        return (-1);
   } else {
      input_filename = argv[1];
      if (argc == 4) {
         restart_read_dirname = argv[2];
         restore_num = atoi(argv[3]);

         is_from_restart = true;
      }
   }

   tbox::plog << "input_filename = " << input_filename << endl;
   tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
   tbox::plog << "restore_num = " << restore_num << endl;

   
   //load input file into an InputDatabase => input_db
   tbox::Pointer<tbox::Database> input_db(new tbox::InputDatabase("input_db"));
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   //grab the global inputs section of the input database and set values
   if (input_db->keyExists("GlobalInputs")) {
      tbox::Pointer<tbox::Database> global_db = 
         input_db->getDatabase("GlobalInputs");
#ifdef SGS
      if (global_db->keyExists("tag_clustering_method")) {
         string tag_clustering_method = 
            global_db->getString("tag_clustering_method");
         mesh::BergerRigoutsos::setClusteringOption(tag_clustering_method); 
      } 
#endif
      // if (global_db->keyExists("refine_schedule_generation_method")) {
      //    string refine_schedule_generation_method =
      //       global_db->getString("refine_schedule_generation_method");
      //    xfer::RefineSchedule::setScheduleGenerationMethod(
      //                                refine_schedule_generation_method);
      // }
      if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
         bool flag = global_db->
            getBool("call_abort_in_serial_instead_of_exit");
         tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
      }
   }

   //read main section into main_db
   tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

	const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));


   const std::string base_name =
      main_db->getStringWithDefault("base_name", "unnamed");

   const std::string log_filename =
      main_db->getStringWithDefault("log_filename", base_name + ".log");

   bool log_all_nodes = false;
   if (main_db->keyExists("log_all_nodes")) {
      log_all_nodes = main_db->getBool("log_all_nodes");
   }
   if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_filename);
   } else {
      tbox::PIO::logOnlyNodeZero(log_filename);
   }

   //viz type and interval
   int viz_dump_interval = 0;
   if (main_db->keyExists("viz_dump_interval")) {
      viz_dump_interval = main_db->getInteger("viz_dump_interval");
   }

   const std::string visit_dump_dirname =
      main_db->getStringWithDefault("viz_dump_dirname", base_name + ".visit");

   int visit_number_procs_per_file = 1;
   if (viz_dump_interval > 0) {
      if (main_db->keyExists("visit_number_procs_per_file")) {
         visit_number_procs_per_file =
            main_db->getInteger("visit_number_procs_per_file");
      }
   }

   //matlab pencil index and interval and filename
   string matlab_dump_filename;
   string matlab_dump_dirname;
   int matlab_dump_interval = 0;
   int  matlab_pencil_direction = 0;
   tbox::Array<int> matlab_pencil_index(dim.getValue()-1);
   bool matlab_default_pencil = true; 
   for (int id = 0; id < dim.getValue()-1; id++) {
      matlab_pencil_index[id] = 0;
   }
   
   if (main_db->keyExists("matlab_dump_interval")) {
      matlab_dump_interval = main_db->getInteger("matlab_dump_interval");
   }
   if (matlab_dump_interval > 0) {
      if (main_db->keyExists("matlab_dump_filename")) {
         matlab_dump_filename = main_db->getString("matlab_dump_filename");
      }
      if (main_db->keyExists("matlab_dump_dirname")) {
         matlab_dump_dirname = main_db->getString("matlab_dump_dirname");
      }
      if (main_db->keyExists("matlab_pencil_direction")) {
         matlab_pencil_direction = 
            main_db->getInteger("matlab_pencil_direction");
      }
      if (main_db->keyExists("matlab_pencil_index")) {
         matlab_default_pencil = false;
         matlab_pencil_index = main_db->getIntegerArray("matlab_pencil_index");
         if (matlab_pencil_index.getSize() != dim.getValue()-1) {
            TBOX_ERROR("`matlab_pencil_index' has "
                       << matlab_pencil_index.getSize() << " values in input. "
                       << dim.getValue()-1 << " values must be specified when default"
                       << " is overridden.");
         }
      }
   }

   string probe_dump_filename;
   string probe_dump_dirname;
   int probe_dump_interval=0;
   double probe_loc_x=0;
   double probe_loc_y=0;
   if (main_db->keyExists("probe_dump_interval")) {
      probe_dump_interval = main_db->getInteger("probe_dump_interval");
   }
   if (probe_dump_interval > 0) {
      if (main_db->keyExists("probe_dump_filename")) {
         probe_dump_filename = main_db->getString("probe_dump_filename");
      }
      if (main_db->keyExists("probe_dump_dirname")) {
         probe_dump_dirname = main_db->getString("probe_dump_dirname");
      }
      if (main_db->keyExists("probe_loc_x")) {
		 probe_loc_x = main_db->getDouble("probe_loc_x");
      }
      if (main_db->keyExists("probe_loc_y")) {
		 probe_loc_y = main_db->getDouble("probe_loc_y");
      }
   }

   //restart interval and directory
   int restart_interval = 0;
   if (main_db->keyExists("restart_interval")) {
      restart_interval = main_db->getInteger("restart_interval");
   }

   const std::string restart_write_dirname =
      main_db->getStringWithDefault("restart_write_dirname",
         base_name + ".restart");

   bool use_refined_timestepping = true;
   if (main_db->keyExists("timestepping")) {
      string timestepping_method = main_db->getString("timestepping");
      if (timestepping_method == "SYNCHRONIZED") {
         use_refined_timestepping = false;
      }
   }

#if (TESTING == 1) && !(HAVE_HDF5)
   /* 
    * If we are autotesting on a system w/o HDF5, the read from
    * restart will result in an error.  We want this to happen
    * for users, so they know there is a problem with the restart,
    * but we don't want it to happen when autotesting.
    */
   is_from_restart  = false;
   restart_interval = 0;
#endif

   //if we will be writing restart files, start the manager
       const bool write_restart = (restart_interval > 0)
         && !(restart_write_dirname.empty());

      /*
       * Get restart manager and root restart database.  If run is from
       * restart, open the restart file.
       */

      tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

#ifdef HAVE_SILO
      /*
       * If SILO is present then use SILO as the file storage format
       * for this example, otherwise it will default to HDF5.
       */
      tbox::Pointer<tbox::SiloDatabaseFactory> silo_database_factory(
         new tbox::SiloDatabaseFactory());
      restart_manager->setDatabaseFactory(silo_database_factory);
#endif

      if (is_from_restart) {
         restart_manager->
         openRestartFile(restart_read_dirname, restore_num,
            mpi.getSize());
      }

   //start the timer manager, timers are listed in the input file
   tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));


   //-------------------------------------------------------------
   // create major algorithm and data objects
   //-------------------------------------------------------------

   //Cartesian grid geometry => grid_geometry
   tbox::Pointer<geom::CartesianGridGeometry> grid_geometry(
       new geom::CartesianGridGeometry(dim,
          "CartesianGeometry",
          input_db->getDatabase("CartesianGeometry")));

   //Patch hierarchy => patch_hierarchy
     tbox::Pointer<hier::PatchHierarchy> patch_hierarchy(
         new hier::PatchHierarchy("PatchHierarchy", grid_geometry,
            input_db->getDatabase("PatchHierarchy")));

   //the SWE class instance (constructor defined in swe.C)  => swe_model
   swe* swe_model = new swe("swe",dim,
                                  input_db->getDatabase("swe"),
                                  grid_geometry);
  
   //Hyperbolic level integrator  => hyp_level_integrator
	tbox::Pointer<algs::HyperbolicLevelIntegrator> hyp_level_integrator(
		new algs::HyperbolicLevelIntegrator("HyperbolicLevelIntegrator",
		input_db->getDatabase(
		"HyperbolicLevelIntegrator"),
		swe_model, true,
		use_refined_timestepping));

   //error detection => error_detector
    tbox::Pointer<mesh::StandardTagAndInitialize> error_detector(
         new mesh::StandardTagAndInitialize(dim,
			"StandardTagAndInitialize",
			hyp_level_integrator,
			input_db->getDatabase("StandardTagAndInitialize")));

   //box generator => box_generator
  		tbox::Pointer<mesh::BergerRigoutsos> box_generator(
         new mesh::BergerRigoutsos(
			dim,
			input_db->getDatabaseWithDefault(
			"BergerRigoutsos",
 			SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL))));

   //load balancer => load_balancer
		tbox::Pointer<mesh::TreeLoadBalancer> load_balancer(
			new mesh::TreeLoadBalancer(dim,
			"LoadBalancer", input_db->getDatabase("LoadBalancer")));
			load_balancer->setSAMRAI_MPI(
         SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

   //gridding => gridding_algorithm
      tbox::Pointer<mesh::GriddingAlgorithm> gridding_algorithm(
         new mesh::GriddingAlgorithm(
			patch_hierarchy,
			"GriddingAlgorithm",
			input_db->getDatabase("GriddingAlgorithm"),
			error_detector,
			box_generator,
			load_balancer));

   //time refinement => time_refinement 
      tbox::Pointer<algs::TimeRefinementIntegrator> time_integrator(
         new algs::TimeRefinementIntegrator("TimeRefinementIntegrator",
			input_db->getDatabase(
			"TimeRefinementIntegrator"),
			patch_hierarchy,
			hyp_level_integrator,
 			gridding_algorithm));

   //visualization 
#ifdef HAVE_HDF5
      tbox::Pointer<appu::VisItDataWriter> visit_data_writer(
         new appu::VisItDataWriter(dim,
            "SWE VisIt Writer",
            visit_dump_dirname,
            visit_number_procs_per_file));
      swe_model->registerVisItDataWriter(visit_data_writer);
#endif

  
   //initialize hierarchy
   double dt_now = time_integrator->initializeHierarchy();

   tbox::RestartManager::getManager()->closeRestartFile();


#if (TESTING == 1)
   /*
    * Create the autotesting component which will verify correctness
    * of the problem. If no automated testing is done, the object does 
    * not get used.
    */
   AutoTester autotester("AutoTester", input_db);
#endif


   //dump SWE class data to log (printClassData defined in swe.C)
#if 1
	tbox::plog << "\nCheck input data and variables before simulation:"
                 << endl;
	tbox::plog << "Input database..." << endl;
	input_db->printClassData(tbox::plog);
	tbox::plog << "\nVariable database..." << endl;
	hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

#endif
	tbox::plog << "\nCheck SWE data... " << endl;
	swe_model->printClassData(tbox::plog);

   //timers to measure I/O
   tbox::Pointer<tbox::Timer> t_write_viz = tbox::TimerManager::getManager()->
                                getTimer("apps::main::write_viz");
   tbox::Pointer<tbox::Timer> t_write_restart = tbox::TimerManager::getManager()->
                                    getTimer("apps::main::write_restart");

   //write initial condition file
		t_write_viz->start();
      if (matlab_dump_interval > 0) {
         dumpMatlabData1dPencil(matlab_dump_dirname,
            matlab_dump_filename,
            0,
            time_integrator->getIntegratorTime(),
            patch_hierarchy,
            matlab_pencil_direction,
            matlab_default_pencil,
            matlab_pencil_index,
            swe_model);
      }
#ifdef HAVE_HDF5
      if (viz_dump_interval > 0) {

         visit_data_writer->writePlotData(
            patch_hierarchy,
            time_integrator->getIntegratorStep(),
            time_integrator->getIntegratorTime());
      }
#endif
      t_write_viz->stop();


#ifdef OPENCL
	//initialize OpenCL environment
  cl_platform_id cpPlatform;
  clGetPlatformIDs(1, &cpPlatform, NULL);

  // Get a GPU device
  cl_device_id cdDevice;
  clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, 1, &cdDevice, NULL);

  char cBuffer[1024];
  clGetDeviceInfo(cdDevice, CL_DEVICE_NAME, sizeof(cBuffer), &cBuffer, NULL);
  printf("CL_DEVICE_NAME:       %s\n", cBuffer);
  clGetDeviceInfo(cdDevice, CL_DRIVER_VERSION, sizeof(cBuffer), &cBuffer, NULL);
  printf("CL_DRIVER_VERSION: %s\n\n", cBuffer);
#endif


   //timing 
   double loop_time = time_integrator->getIntegratorTime();
   double loop_time_end = time_integrator->getEndTime();

   // initialize matlab frame counter
   int matlab_cnt = 1;

   // initialize probe frame counter
   int probe_cnt = 1;

   //-------------------------------------------------------------
   // BEGIN: main integration loop over time
   //-------------------------------------------------------------
   while ( (loop_time < loop_time_end) && 
          time_integrator->stepsRemaining() ) {

      int iteration_num = time_integrator->getIntegratorStep() + 1;

      //report progress to logfile
      tbox::plog <<endl<<endl;
      //tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
	  tbox::pout << "Timestep " << iteration_num   << " Time " ;
	  tbox:: pout << loop_time << " DT " << dt_now << endl;
      //tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      //tbox::plog <<endl<<endl;

      //integrate
      double dt_new = time_integrator->advanceHierarchy(dt_now);

      //update time
      loop_time += dt_now;
      dt_now = dt_new;

      //report progress
      // tbox::plog <<endl<<endl;
      //      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      //      tbox::pout << "At end of timestep # " <<  iteration_num - 1 << endl;
      //      tbox::pout << "Simulation time is " << loop_time << endl;
      //      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      //      tbox::plog <<endl<<endl;

      //dump restart
      if ( write_restart ) {

         if ( (iteration_num % restart_interval) == 0 ) { 
            t_write_restart->start();
            tbox::RestartManager::getManager()->
               writeRestartFile(restart_write_dirname,
                                iteration_num);
            t_write_restart->stop();
         }
      }

      //dump output
      t_write_viz->start();
      if ( (viz_dump_interval > 0) 
           && (iteration_num % viz_dump_interval) == 0 ) {
#ifdef HAVE_HDF5
         if (uses_visit) {
            visit_data_writer->writePlotData(patch_hierarchy,
                                             iteration_num,
                                             loop_time);
         }
#endif
         
      }
      if ((matlab_dump_interval > 0)
           && (iteration_num % matlab_dump_interval) == 0 ) {
          dumpMatlabData1dPencil(matlab_dump_dirname,
                                 matlab_dump_filename, 
                                 matlab_cnt,
                                 loop_time,
                                 patch_hierarchy,
                                 matlab_pencil_direction,
                                 matlab_default_pencil,
                                 matlab_pencil_index,
                                 swe_model);
		 matlab_cnt++;
      }
      if ((probe_dump_interval > 0)
           && (iteration_num % probe_dump_interval) == 0 ) {
          dumpProbe(probe_dump_dirname,
                                 probe_dump_filename, 
                                 probe_cnt,
                                 loop_time,
                                 patch_hierarchy,
                                 probe_loc_x,
                                 probe_loc_x,
                                 swe_model);
		 probe_cnt++;
      }
      t_write_viz->stop();


   } 
   //-------------------------------------------------------------
   // END: main integration loop over time
   //-------------------------------------------------------------

   //dump timer info
   tbox::TimerManager::getManager()->print(tbox::plog);

   //nullify all objects
   patch_hierarchy.setNull();
   grid_geometry.setNull();

   box_generator.setNull();
   load_balancer.setNull();
   hyp_level_integrator.setNull();
   error_detector.setNull();
   gridding_algorithm.setNull();
   time_integrator.setNull();
#ifdef HAVE_HDF5
   visit_data_writer.setNull();
#endif

   if (swe_model) delete swe_model;

   }

   if (num_failures == 0) {
      tbox::pout << "\nPASSED:  swe" << endl;
   }

   //shutdown
   tbox::SAMRAIManager::shutdown();
	tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return(num_failures); 
}

static void dumpMatlabData1dPencil(
   const string& dirname,
   const string& filename,
   const int ext,
   const double plot_time,
   const tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const int pencil_direction,
   const bool default_pencil,
   const tbox::Array<int>& pencil_index,
   swe* swe_model)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   /*
    * Compute the boxes to write out data at each level of the hierarchy.
    */

   const tbox::Dimension& dim = hierarchy->getDim();

   int nlevels = 1;

   if ((dim == tbox::Dimension(1))) {
      nlevels = hierarchy->getNumberOfLevels();
   }

   hier::BoxList domain(hierarchy->getGridGeometry()->getPhysicalDomain());
   hier::Box pencil_box(domain.getBoundingBox());

   if (dim > tbox::Dimension(1)) {
      int indx = 0;
      int id = 0;
      tbox::Array<int> tmp(dim.getValue() - 1);
      for (id = 0; id < dim.getValue() - 1; id++) {
         tmp[id] = pencil_index[id];
      }
      if (default_pencil) {
         hier::Index ifirst = domain.getBoundingBox().lower();
         indx = 0;
         for (id = 0; id < dim.getValue(); id++) {
            if (id != pencil_direction) {
               tmp[indx] = ifirst(id);
               indx++;
            }
         }
      }
      indx = 0;
      for (id = 0; id < dim.getValue(); id++) {
         if (id != pencil_direction) {
            pencil_box.lower(id) = tmp[indx];
            pencil_box.upper(id) = tmp[indx];
            indx++;
         }
      }
   }

   tbox::Array<hier::BoxList> outboxes(nlevels);

   for (int l1 = 0; l1 < nlevels; l1++) {
      tbox::Pointer<hier::PatchLevel> level = hierarchy->getPatchLevel(l1);
      outboxes[l1] = hier::BoxList(level->getBoxes());

      if (l1 < nlevels - 1) {

         tbox::Pointer<hier::PatchLevel> finer_level =
            hierarchy->getPatchLevel(l1 + 1);
         hier::IntVector coarsen_ratio =
            finer_level->getRatioToCoarserLevel();
         hier::BoxList takeaway = hier::BoxList(finer_level->getBoxes());
         takeaway.coarsen(coarsen_ratio);
         outboxes[l1].removeIntersections(takeaway);
      }

   }

   /*
    * Create matlab filename and open the output stream.
    */

   string dump_filename = filename;
   if (!dirname.empty()) {
      tbox::Utilities::recursiveMkdir(dirname);
      dump_filename = dirname;
      dump_filename += "/";
      dump_filename += filename;
   }

   const int size = static_cast<int>(dump_filename.length()) + 16;
   char* buffer = new char[size];

   if (mpi.getSize() > 1) {
      sprintf(buffer, "%s.%04d.dat.%05d", dump_filename.c_str(),
         ext, mpi.getRank());
   } else {
      sprintf(buffer, "%s_%04d.dat", dump_filename.c_str(), ext);
   }

   /*
    * Open a new output file having name name of buffer character array.
    */

   ofstream outfile(buffer, ios::out);
   outfile.setf(ios::scientific);
   outfile.precision(10);

   delete[] buffer;

   /*
    * There are 7 values dumped for every cell.  Here we dump the time.
    */
   for (int i = 0; i < 6 + 1; i++) {
      outfile << plot_time << "  ";
   }
   outfile << endl;

   swe_model->setDataContext(
      hier::VariableDatabase::getDatabase()->getContext("CURRENT"));

   for (int l5 = 0; l5 < nlevels; l5++) {
      tbox::Pointer<hier::PatchLevel> level = hierarchy->getPatchLevel(l5);

      hier::Box level_pencil_box = pencil_box;
      if (l5 > 0) {
         level_pencil_box.refine(level->getRatioToLevelZero());
      }

      for (hier::PatchLevel::Iterator i(level); i; i++) {
         tbox::Pointer<hier::Patch> patch = *i;
         hier::Box pbox = patch->getBox();

         for (hier::BoxList::Iterator b(outboxes[l5]); b; b++) {
            const hier::Box box = b() * pbox * level_pencil_box;

            swe_model->writeData1dPencil(patch,
               box,
               pencil_direction,
               outfile);
         }

      }

   }

   swe_model->clearDataContext();

   outfile.close();

}


// dumpprobe - determine the finest patch level containing the probe point
// and determing the nearest cell
// call 
void dumpProbe(const string& dirname,
               const string& filename,
               const int ext,
               const double plot_time,
               const tbox::Pointer<hier::PatchHierarchy > hierarchy,
               const double probe_x,
               const double probe_y,
               swe* swe_model)
{
	cout << " somthing\n";
//    /*
//     * Compute the boxes to write out data at each level of the hierarchy.
//     */
// 
//    int nlevels = 1;
// 
// 
// 
//    hier::BoxList<NDIM> domain(hierarchy->getGridGeometry()->getPhysicalDomain());
//    hier::Box<NDIM> probe_box(domain.getBoundingBox());
// 
//    int indx = 0;
//    int id = 0;
//    tbox::Array<int> tmp(NDIM-1);
//    for (id = 0; id < NDIM-1; id++) {
//       tmp[id] = pencil_index[id];
//    }
//    if (default_pencil) {
//       hier::Index<NDIM> ifirst = domain.getBoundingBox().lower();
//       indx = 0;
//       for (id = 0; id < NDIM; id++) {
//          if (id != pencil_direction) { 
//             tmp[indx] = ifirst(id);
//             indx++;
//          }
//       } 
//    }
//    indx = 0;
//    for (id = 0; id < NDIM; id++) {
//       if (id != pencil_direction) {
//          pencil_box.lower(id) = tmp[indx];
//          pencil_box.upper(id) = tmp[indx];
//          indx++;
//       }
//    }
// 
//    tbox::Array<hier::BoxList<NDIM> > outboxes(nlevels);
// 
//    for (int l1 = 0; l1 < nlevels; l1++) {
//       tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(l1);
//       outboxes[l1] = level->getBoxes();
// 
//       if (l1 < nlevels-1) {
// 
// 	 tbox::Pointer<hier::PatchLevel<NDIM> > finer_level = 
// 	    hierarchy->getPatchLevel(l1+1);
//          hier::IntVector<NDIM> coarsen_ratio 
// 	    = finer_level->getRatioToCoarserLevel();
//          hier::BoxList<NDIM> takeaway = finer_level->getBoxes();
//          takeaway.coarsen(coarsen_ratio);
//          outboxes[l1].removeIntersections(takeaway);
//       }
// 
//    }
// 
//    //
//    // Create probe filename and open the output stream.
//    //
// 
//    string dump_filename = filename;
//    if (!dirname.empty()) {
//       tbox::Utilities::recursiveMkdir(dirname);
//       dump_filename = dirname;
//       dump_filename += "/";
//       dump_filename += filename;
//    }
// 
//    const int size = dump_filename.length() + 16;
//    char *buffer = new char[size];
// 
//    if (tbox::SAMRAI_MPI::getNodes() > 1) {
//       sprintf(buffer, "%s.%04d.dat.%05d", dump_filename.c_str(),
//               ext, tbox::SAMRAI_MPI::getRank());
//    }
//    else {
//       sprintf(buffer, "%s_%04d.dat", dump_filename.c_str(), ext);
//    }
// 
//    /*
//     * Open a new output file having name name of buffer character array.
//     */
// 
//    ofstream outfile(buffer, ios::out);
//    outfile.setf(ios::scientific);
//    outfile.precision(10);
// 
//    delete [] buffer;
// 
//    /*
//     * There are 7 values dumped for every cell.  Here we dump the time.
//     */
// //   for (int i = 0; i < 6+1; i++) {
//       outfile << plot_time << "  ";
// //   }
//    outfile << endl;
// 
//    swe_model->setDataContext(
//       hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT") );
// 
//    for (int l5 = 0; l5 < nlevels; l5++) {
//       tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(l5);
// 
//       hier::Box<NDIM> level_pencil_box = pencil_box;
//       if (l5 > 0) {
//          level_pencil_box.refine(level->getRatio());
//       }
// 
//       for (hier::PatchLevel<NDIM>::Iterator i(level); i; i++) {
//          tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(i());
//          hier::Box<NDIM> pbox = patch->getBox();
// 
//          for (hier::BoxList<NDIM>::Iterator b(outboxes[l5]); b; b++) {
//             const hier::Box<NDIM> box = b() * pbox * level_pencil_box; 
// 
// 			//swe_model->writeProbe(patch,ival,jval,outfile);
//                                   
//          }
// 
//       }

//   }

   //close up!
   //swe_model->clearDataContext();

   //outfile.close();

}

