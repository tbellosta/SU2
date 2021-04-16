/*!
 * \file CPTDriver.cpp
 * \brief The main subroutines for driving single-zone problems.
 * \author T. Bellosta
 * \version 7.0.7 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/drivers/CPTDriver.hpp"
#include "../../include/solvers/CPTSolver.hpp"
#include "../../include/definition_structure.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIteration.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

CPTDriver::CPTDriver(char* confFile,
                       unsigned short val_nZone,
                       SU2_Comm MPICommunicator) : CDriver(confFile,
                                                          val_nZone,
                                                          MPICommunicator,
                                                          false) {

  /*--- Initialize the counter for TimeIter ---*/
  TimeIter = 0;

  StopCalc_PT = false;
  runtimeFlow = true;
  runtimeSplashing = false;


}


CPTDriver::CPTDriver(char* confFile,
                       unsigned short val_nZone,
                       SU2_Comm MPICommunicator, bool splashingConfig) : CDriver(confFile,
                                                          val_nZone,
                                                          MPICommunicator,
                                                          false) {

  /*--- Initialize the counter for TimeIter ---*/
  TimeIter = 0;


  StopCalc_PT = false;
  runtimeFlow = true;
  runtimeSplashing = false;
  splashing = splashingConfig;


}


CPTDriver::~CPTDriver(void) {

}

void CPTDriver::StartSolver() {

  StartTime = SU2_MPI::Wtime();

  config_container[ZONE_0]->Set_StartTime(StartTime);

  /*--- Main external loop of the solver. Runs for the number of time steps required. ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;

  if (rank == MASTER_NODE){
    cout << endl <<"Simulation Run using the particle tracking (Single-zone) Driver" << endl;
    if (driver_config->GetTime_Domain())
      cout << "The simulation will run for "
           << driver_config->GetnTime_Iter() - config_container[ZONE_0]->GetRestart_Iter() << " time steps." << endl;
  }

  /*--- Set the initial time iteration to the restart iteration. ---*/
  if (config_container[ZONE_0]->GetRestart() && driver_config->GetTime_Domain())
    TimeIter = config_container[ZONE_0]->GetRestart_Iter();

  /*--- Run the problem until the number of time iterations required is reached. ---*/
  while ( TimeIter < config_container[ZONE_0]->GetnTime_Iter() ) {

    /*--- Solve flow system then particles ---*/

    /*--- Flow system ---*/
    setRuntimeFlow();
    /*--- Perform some preprocessing before starting the time-step simulation. ---*/

    Preprocess(TimeIter);

    /*--- Run a time-step iteration of the single-zone problem. ---*/

    Run();

    /*--- Perform some postprocessing on the solution before the update ---*/

    Postprocess();

    /*--- Update the solution for dual time stepping strategy ---*/

    Update();

    /*--- Monitor the computations after each iteration. ---*/

    Monitor(TimeIter);

    /*--- Output the solution in files. ---*/

    Output(TimeIter);

    /*--- Particle system ---*/
    setRuntimeParticles();


    unsigned short FinestMesh = config_container[ZONE_0]->GetFinestMesh(); 
    CSolver** solvers_fine = solver_container[ZONE_0][INST_0][FinestMesh];
    CPTSolver *PTsolver = dynamic_cast<CPTSolver*>(solvers_fine[PT_SOL]);
    su2double* MVD_multibin = config_container[ZONE_0]->GetMVDMultibin();
    su2double* percentage_multibin = config_container[ZONE_0]->GetPercentageMultibin();
    su2double* CFL_multibin = config_container[ZONE_0]->GetCFLMultibin();
    unsigned short nBins = config_container[ZONE_0]->GetNBins();
    bool multibin = config_container[ZONE_0]-> GetMultiBin();
    PTsolver->multibin=multibin;
    su2double LWC_inf = config_container[ZONE_0]->GetLiquidWaterContent();
    


    if(multibin){ //only steady multibin

      int i = 0;  

      if(rank==MASTER_NODE){
        cout << "\n\n RUNNING PARTICLE TRACKING WITH MULTIBIN DATA: \n";
        for (i = 0;i<nBins;i++){
          cout << " "<<i+1<<")\t " << percentage_multibin[i]<< "%   \t-    MVD = "<< MVD_multibin[i] << " micrometers\n";
        
        }
      }

      PTsolver->InitializeMultiBin(MVD_multibin, percentage_multibin, CFL_multibin, LWC_inf, nBins);

      //running simulation for each bin
      for (i = 0;i<nBins;i++){
        if(rank==MASTER_NODE){
          cout << "\n\nRUNNING BIN "<<i+1<<")\t" << percentage_multibin[i]<< "%   \t-    MVD = "<< MVD_multibin[i] << " micrometers\n\n";
        }

        //not sure i should simulate each bin with the same LWC_inf
        PTsolver->SetBin(i);
        
        Preprocess(TimeIter);
        Run();
        Postprocess();
        Update();
        Monitor(TimeIter);
        #ifdef HAVE_MPI
          SU2_MPI::Barrier(MPI_COMM_WORLD);
        #endif
        if(splashing){   
            CGeometry* geometry_fine = geometry_container[ZONE_0][INST_0][FinestMesh];
            //CSolver** solvers_fine = solver_container[ZONE_0][INST_0][FinestMesh];
            CPTSolver *splashingsolver = dynamic_cast<CPTSolver*>(solvers_fine[SPLASHINGPT_SOL]);
            //CPTSolver *PTsolver = dynamic_cast<CPTSolver*>(solvers_fine[PT_SOL]);

            /*--- Sum this bin's collection efficiency ---*/
            PTsolver->AddBinCollectionEfficiency(geometry_fine);
            /*--- Compute BCs for splashing droplets (and corrects BETA) ---*/
            PTsolver->ComputeSplashingBCs(geometry_fine,splashingsolver,config_container[ZONE_0],runtimeSplashing); 
        }
        #ifdef HAVE_MPI
          SU2_MPI::Barrier(MPI_COMM_WORLD);
        #endif
        Output(TimeIter);

        #ifdef HAVE_MPI
          SU2_MPI::Barrier(MPI_COMM_WORLD);
        #endif
      }
    }
    else{

      /*--- Perform some preprocessing before starting the time-step simulation. ---*/

      Preprocess(TimeIter);

      /*--- Run a time-step iteration of the single-zone problem. ---*/

      Run();



      /*--- Perform some postprocessing on the solution before the update ---*/

      Postprocess();

      /*--- Update the solution for dual time stepping strategy ---*/

      Update();

      /*--- Monitor the computations after each iteration. ---*/

      Monitor(TimeIter);

      #ifdef HAVE_MPI
        SU2_MPI::Barrier(MPI_COMM_WORLD);
      #endif

      //cout<<"\n\n"<< rank <<" before splash \n\n";
      if(splashing){   
        CGeometry* geometry_fine = geometry_container[ZONE_0][INST_0][FinestMesh];
        //CSolver** solvers_fine = solver_container[ZONE_0][INST_0][FinestMesh];
        CPTSolver *splashingsolver = dynamic_cast<CPTSolver*>(solvers_fine[SPLASHINGPT_SOL]);
        //CPTSolver *PTsolver = dynamic_cast<CPTSolver*>(solvers_fine[PT_SOL]);

        
        /*--- Compute BCs for splashing droplets (and corrects BETA) ---*/
        //cout<<"\n\n\t"<< rank <<" before ComputeSplashingBCs \n\n";
        PTsolver->ComputeSplashingBCs(geometry_fine,splashingsolver,config_container[ZONE_0],runtimeSplashing);
          
        //cout<<"\n\n\t"<< rank <<" after ComputeSplashingBCs \n\n";  
      }
      //cout<<"\n\n"<< rank <<" before driver->output() \n\n";

      #ifdef HAVE_MPI
        SU2_MPI::Barrier(MPI_COMM_WORLD);
      #endif
      /*--- Output the solution in files. ---*/
      Output(TimeIter);

      #ifdef HAVE_MPI
        SU2_MPI::Barrier(MPI_COMM_WORLD);
      #endif

  //     if(splashing){
  //       /*--- Set splashing as runtime ---*/
  //       setRuntimeSplashing();

  //       /*--- Compute BCs for splashing droplets ---*/
        
  //       unsigned short FinestMesh = config_container[ZONE_0]->GetFinestMesh();    
  //       CGeometry* geometry_fine = geometry_container[ZONE_0][INST_0][FinestMesh];
  //       CSolver** solvers_fine = solver_container[ZONE_0][INST_0][FinestMesh];
  //       CPTSolver *splashingsolver = dynamic_cast<CPTSolver*>(solvers_fine[SPLASHINGPT_SOL]);
  //       CPTSolver *PTsolver = dynamic_cast<CPTSolver*>(solvers_fine[PT_SOL]);

        
  //       /*--- Compute BCs for splashing droplets (and corrects BETA) ---*/
  //       PTsolver->ComputeSplashingBCs(geometry_fine,splashingsolver,config_container[ZONE_0],runtimeSplashing);
            
        


  // #ifdef HAVE_MPI
  //   SU2_MPI::Barrier(MPI_COMM_WORLD);
  // #endif

  //       cout << "\n before preprocess \n";
  //       Preprocess(TimeIter);

  //       /*--- Run a time-step iteration of the single-zone problem. ---*/
  //       cout << "\n before run \n";
  //       Run();

  //       /*--- Perform some postprocessing on the solution before the update ---*/
  //       cout << "\n before postprocess \n";
  //       Postprocess();

  //       /*--- Update the solution for dual time stepping strategy ---*/

  //       Update();

  //       /*--- Monitor the computations after each iteration. ---*/

  //       Monitor(TimeIter);

  //       /*--- Output the solution in files. ---*/
  //       Output(TimeIter);

  //     }

      
      


      /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    }

    if (StopCalc) break;

    TimeIter++;

  }

}


void CPTDriver::Preprocess(unsigned long TimeIter) {
  



  /*--- Set runtime option ---*/


  Runtime_Options();

  /*--- Set the current time iteration in the config ---*/

  config_container[ZONE_0]->SetTimeIter(TimeIter);

  /*--- Store the current physical time in the config container, as
   this can be used for verification / MMS. This should also be more
   general once the drivers are more stable. ---*/

  if (config_container[ZONE_0]->GetTime_Marching())
    config_container[ZONE_0]->SetPhysicalTime(static_cast<su2double>(TimeIter)*config_container[ZONE_0]->GetDelta_UnstTimeND());
  else
    config_container[ZONE_0]->SetPhysicalTime(0.0);

  /*--- Set the initial condition for EULER/N-S/RANS ---------------------------------------------*/
  if (runtimeFlow) {
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[ZONE_0][INST_0],
                                                                            solver_container[ZONE_0][INST_0],
                                                                            config_container[ZONE_0], TimeIter);
  }
  else {
    if(runtimeSplashing) {

      solver_container[ZONE_0][INST_0][MESH_0][SPLASHINGPT_SOL]->SetInitialCondition(
      geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], TimeIter);
      
    } 
    else {
      solver_container[ZONE_0][INST_0][MESH_0][PT_SOL]->SetInitialCondition(
      geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], TimeIter);
    }
  }


#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif


  /*--- Run a predictor step ---*/
  if (config_container[ZONE_0]->GetPredictor() && runtimeFlow)
    iteration_container[ZONE_0][INST_0]->Predictor(output_container[ZONE_0], integration_container, geometry_container, solver_container,
        numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Perform a dynamic mesh update if required. ---*/
  /*--- For the Disc.Adj. of a case with (rigidly) moving grid, the appropriate
          mesh cordinates are read from the restart files. ---*/
  if (!(config_container[ZONE_0]->GetGrid_Movement() && config_container[ZONE_0]->GetDiscrete_Adjoint()))
    DynamicMeshUpdate(TimeIter);

}

void CPTDriver::Run() {

  unsigned long OuterIter = 0;
  config_container[ZONE_0]->SetOuterIter(OuterIter);

  /*--- Iterate the zone as a block, either to convergence or to a max number of iterations ---*/
  if (runtimeFlow) {
    iteration_container[ZONE_0][INST_0]->Solve(output_container[ZONE_0], integration_container, geometry_container,
                                               solver_container, numerics_container, config_container, surface_movement,
                                               grid_movement, FFDBox, ZONE_0, INST_0);
  } else {
    if(runtimeSplashing){

      iteration_container_splashingPT[ZONE_0][INST_0]->Solve(output_container_splashingPT[ZONE_0], integration_container, geometry_container,
                                               solver_container, numerics_container, config_container, surface_movement,
                                               grid_movement, FFDBox, ZONE_0, INST_0);  

    }
    else{
      iteration_container_PT[ZONE_0][INST_0]->Solve(output_container_PT[ZONE_0], integration_container, geometry_container,
                                               solver_container, numerics_container, config_container, surface_movement,
                                               grid_movement, FFDBox, ZONE_0, INST_0);  

    }
  }

}

void CPTDriver::Postprocess() {

  if (runtimeFlow) {
    iteration_container[ZONE_0][INST_0]->Postprocess(
        output_container[ZONE_0], integration_container, geometry_container, solver_container, numerics_container,
        config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);
  } else {
    if(runtimeSplashing){

      iteration_container_PT[ZONE_0][INST_0]->Postprocess(
        output_container_splashingPT[ZONE_0], integration_container, geometry_container, solver_container, numerics_container,
        config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    }
    else{

      iteration_container_PT[ZONE_0][INST_0]->Postprocess(
        output_container_PT[ZONE_0], integration_container, geometry_container, solver_container, numerics_container,
        config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    }
  }

    /*--- A corrector step can help preventing numerical instabilities ---*/

    if (config_container[ZONE_0]->GetRelaxation() && runtimeFlow)
      iteration_container[ZONE_0][INST_0]->Relaxation(output_container[ZONE_0], integration_container, geometry_container, solver_container,
          numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

}

void CPTDriver::Update() {

  if (runtimeFlow) {
    iteration_container[ZONE_0][INST_0]->Update(output_container[ZONE_0], integration_container, geometry_container,
                                                solver_container, numerics_container, config_container,
                                                surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);
  } else {
    if(runtimeSplashing){
      iteration_container_PT[ZONE_0][INST_0]->Update(output_container_splashingPT[ZONE_0], integration_container, geometry_container,
                                                solver_container, numerics_container, config_container,
                                                surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    }else{
      iteration_container_PT[ZONE_0][INST_0]->Update(output_container_PT[ZONE_0], integration_container, geometry_container,
                                                solver_container, numerics_container, config_container,
                                                surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    }
  }

}

void CPTDriver::Output(unsigned long TimeIt) {

  /*--- Time the output for performance benchmarking. ---*/

  
  StopTime = SU2_MPI::Wtime();

  UsedTimeCompute += StopTime-StartTime;

  StartTime = SU2_MPI::Wtime();

  bool wrote_files;
  if (runtimeFlow)
    wrote_files = output_container[ZONE_0]->SetResult_Files(geometry_container[ZONE_0][INST_0][MESH_0],
                                                               config_container[ZONE_0],
                                                               solver_container[ZONE_0][INST_0][MESH_0],
                                                               TimeIt, StopCalc);
  else{
    if(runtimeSplashing){
      wrote_files = output_container_splashingPT[ZONE_0]->SetResult_Files(geometry_container[ZONE_0][INST_0][MESH_0],
                                                               config_container[ZONE_0],
                                                               solver_container[ZONE_0][INST_0][MESH_0],
                                                               TimeIt, StopCalc);

    }else{
      wrote_files = output_container_PT[ZONE_0]->SetResult_Files(geometry_container[ZONE_0][INST_0][MESH_0],
                                                               config_container[ZONE_0],
                                                               solver_container[ZONE_0][INST_0][MESH_0],
                                                               TimeIt, StopCalc);

    }
  }

  if (wrote_files){

    StopTime = SU2_MPI::Wtime();

    UsedTimeOutput += StopTime-StartTime;
    OutputCount++;
    BandwidthSum = config_container[ZONE_0]->GetRestart_Bandwidth_Agg();

    StartTime = SU2_MPI::Wtime();

    config_container[ZONE_0]->Set_StartTime(StartTime);
  }

}

void CPTDriver::DynamicMeshUpdate(unsigned long TimeIter) {

  CIteration *iteration = (runtimeFlow) ? iteration_container[ZONE_0][INST_0] : iteration_container_PT[ZONE_0][INST_0];

  /*--- Legacy dynamic mesh update - Only if GRID_MOVEMENT = YES ---*/
  if (config_container[ZONE_0]->GetGrid_Movement()) {
    iteration->SetGrid_Movement(geometry_container[ZONE_0][INST_0],surface_movement[ZONE_0],
                                grid_movement[ZONE_0][INST_0], solver_container[ZONE_0][INST_0],
                                config_container[ZONE_0], 0, TimeIter);
  }

  /*--- New solver - all the other routines in SetGrid_Movement should be adapted to this one ---*/
  /*--- Works if DEFORM_MESH = YES ---*/
  iteration->SetMesh_Deformation(geometry_container[ZONE_0][INST_0],
                                 solver_container[ZONE_0][INST_0][MESH_0],
                                 numerics_container[ZONE_0][INST_0][MESH_0],
                                 config_container[ZONE_0], NONE);

  /*--- Update the wall distances if the mesh was deformed. ---*/
  if (config_container[ZONE_0]->GetGrid_Movement() ||
      config_container[ZONE_0]->GetDeform_Mesh()) {
    CGeometry::ComputeWallDistance(config_container, geometry_container);
  }
}

bool CPTDriver::Monitor(unsigned long TimeIter){

  unsigned long nInnerIter, InnerIter, nTimeIter;
  su2double MaxTime, CurTime;
  bool TimeDomain, InnerConvergence, TimeConvergence, FinalTimeReached, MaxIterationsReached;

  COutput *output = nullptr;
  if(runtimeFlow){
    output = output_container[ZONE_0];
  }
  else{
    if(runtimeSplashing){
      output = output_container_splashingPT[ZONE_0];
    }
    else{
      output = output_container_PT[ZONE_0];
    }
  }

  nInnerIter = config_container[ZONE_0]->GetnInner_Iter();
  InnerIter  = config_container[ZONE_0]->GetInnerIter();
  nTimeIter  = config_container[ZONE_0]->GetnTime_Iter();
  MaxTime    = config_container[ZONE_0]->GetMax_Time();
  CurTime    = output->GetHistoryFieldValue("CUR_TIME");

  TimeDomain = config_container[ZONE_0]->GetTime_Domain();


  /*--- Check whether the inner solver has converged --- */

  if (TimeDomain == NO){

    InnerConvergence     = output->GetConvergence();
    MaxIterationsReached = InnerIter+1 >= nInnerIter;

    if ((MaxIterationsReached || InnerConvergence) && (rank == MASTER_NODE)) {
      cout << endl << "----------------------------- Solver Exit -------------------------------" << endl;
      if (InnerConvergence) cout << "All convergence criteria satisfied." << endl;
      else cout << endl << "Maximum number of iterations reached (ITER = " << nInnerIter << ") before convergence." << endl;
      output->PrintConvergenceSummary();
      cout << "-------------------------------------------------------------------------" << endl;
    }

    StopCalc = MaxIterationsReached || InnerConvergence;
  }



  if (TimeDomain == YES) {

    /*--- Check whether the outer time integration has reached the final time ---*/

    TimeConvergence = GetTimeConvergence();
    FinalTimeReached     = CurTime >= MaxTime;
    MaxIterationsReached = TimeIter+1 >= nTimeIter;

    if ((FinalTimeReached || MaxIterationsReached || TimeConvergence) && (rank == MASTER_NODE)){
      cout << endl << "----------------------------- Solver Exit -------------------------------";
      if (TimeConvergence)     cout << endl << "All windowed time-averaged convergence criteria are fullfilled." << endl;
      if (FinalTimeReached)     cout << endl << "Maximum time reached (MAX_TIME = " << MaxTime << "s)." << endl;
      if (MaxIterationsReached) cout << endl << "Maximum number of time iterations reached (TIME_ITER = " << nTimeIter << ")." << endl;
      cout << "-------------------------------------------------------------------------" << endl;
    }
    StopCalc = FinalTimeReached || MaxIterationsReached|| TimeConvergence;
  }

  /*--- Reset the inner convergence --- */

  output->SetConvergence(false);

  /*--- Increase the total iteration count --- */

  IterCount += config_container[ZONE_0]->GetInnerIter()+1;

  return StopCalc;
}

void CPTDriver::Runtime_Options(){


  ifstream runtime_configfile;

  /*--- Try to open the runtime config file ---*/

  runtime_configfile.open(runtime_file_name, ios::in);

  /*--- If succeeded create a temporary config object ---*/

  if (runtime_configfile.good()){
    CConfig *runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);
    delete runtime;
  }

}

bool CPTDriver::GetTimeConvergence() const{
  COutput *output = nullptr;
  if(runtimeFlow){
    output = output_container[ZONE_0];

  }
  else{
    if(runtimeSplashing){
      output =  output_container_splashingPT[ZONE_0];
    }
    else{
       output = output_container_PT[ZONE_0];
    }
  }
  return output->GetTimeConvergence();
}
