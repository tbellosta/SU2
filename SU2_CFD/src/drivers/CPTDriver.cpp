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


    int i = 0; 
    unsigned short FinestMesh = config_container[ZONE_0]->GetFinestMesh(); 
    CSolver** solvers_fine = solver_container[ZONE_0][INST_0][FinestMesh];
    CPTSolver *PTsolver = dynamic_cast<CPTSolver*>(solvers_fine[PT_SOL]);
    su2double* diameters_multibin = config_container[ZONE_0]->GetMVDMultibin();
    su2double* percentage_multibin = config_container[ZONE_0]->GetPercentageMultibin();
    su2double* CFL_multibin = config_container[ZONE_0]->GetCFLMultibin();
    unsigned short nBins = 0;

    bool multibin = config_container[ZONE_0]-> GetMultiBin();
    PTsolver->multibin=multibin;

    if(multibin){
      if(diameters_multibin != nullptr && percentage_multibin != nullptr){
        //bins already inserted by user in .cfg
        nBins = config_container[ZONE_0]->GetNBins();
      }
      else{
        //must compute multibin distribution using Langmuir D distribution double 12th order polynomial fit
        nBins = config_container[ZONE_0]->GetNBinsUser();

        if(nBins>60){
          nBins=60;
          if(rank==MASTER_NODE){
            cout <<"+-----------------------------------------------------------------------+" << endl;
            cout<<"|\n|\n| Number of bins requested is too high, maximum is 60. Automatically limited. \n|\n|\n";
            cout <<"-----------------------------------------------------------------------+" << endl;
          }
        }
        else{
          if(nBins<=0){
            nBins=1;
            if(rank==MASTER_NODE){
            cout <<"+-----------------------------------------------------------------------+" << endl;
            cout<<"|\n|\n| Number of bins requested is too low, minimum is 1. Automatically raised. |\n|\n|\n";
            cout <<"+-----------------------------------------------------------------------+" << endl;
            }
          }
        }
        su2double diameters_tmp[nBins];
        su2double percentage_tmp[nBins];

        su2double poly_coeff1[13] = {0.9997,-0.0073,0.8022,-16.1160,129.4137,-621.6733,1805.9536, 
        -3120.2903,2886.3084,-740.7042,-1030.8815,952.6488,-245.9530};
        su2double poly_coeff2[13] = {31.183098532270120,-1.941145818967928e+02,5.442444544336961e+02,-8.664362455492004e+02,
        8.406650710006957e+02,-4.828259304177341e+02,1.187099872116714e+02, 
        41.801778610435974,-48.806507186880935,20.101826063061310,
        -4.555851892681544,0.562876844773370,-0.029765667278845};
        
        //MVD step ratio
        su2double diameter_ratio_step = (su2double)3 / (su2double)nBins;
        su2double MVD = config_container[ZONE_0]->GetParticle_Size();

        su2double previous_LWC_fraction = 1;

        //MVD ratio above which polynomial 1 must be switched for polynomial 2
        su2double diameter_ratio_cutoff12 = 0.9;


        for(i=0; i<nBins;i++){
          su2double diameter_ratio_i = (i+1)* diameter_ratio_step;
          diameters_tmp[i] = diameter_ratio_i * MVD;
          su2double *poly_coeff = (su2double*)malloc(sizeof(su2double)*13);;

          su2double LWC_fraction_i = 0;
          unsigned short n = 0;
          if(diameter_ratio_i < diameter_ratio_cutoff12){
            //use polynomial 1
            unsigned short j=0;
            for(j=0;j<13;j++){
              poly_coeff[j] = poly_coeff1[j];
            }
          }else{
            //use polynomial 2
            unsigned short j=0;
            for(j=0;j<13;j++){
              poly_coeff[j] = poly_coeff2[j];
            }
          }
          //if(rank == MASTER_NODE){unsigned short j=0;
          //  for(j=0;j<13;j++){
          //    cout << "\n poly = " <<poly_coeff[j] << "";
          //  }
          //}


          for(n=0; n<13; n++){
            LWC_fraction_i += poly_coeff[n] * pow(diameter_ratio_i,n);
          }
          if(LWC_fraction_i>=1){
            LWC_fraction_i =1;
          }
          if(LWC_fraction_i<=0){
            LWC_fraction_i = 0;
          }

          //if(rank == MASTER_NODE){
          //  cout << "\n diam_frac = " <<diameter_ratio_i << "";
          //  cout << "\n LWC_frac = " <<LWC_fraction_i << "";
          //  cout << "\n prev_LWC_frac = " <<previous_LWC_fraction << "\n\n";
          //}

          percentage_tmp[i] = 100 * (previous_LWC_fraction - LWC_fraction_i); 
          if(percentage_tmp[i]<0){
            //percentage_tmp[i]<0 = 0;
          }

          

          previous_LWC_fraction = LWC_fraction_i;
        }

        diameters_multibin = (su2double*)malloc(sizeof(su2double)*nBins);
        percentage_multibin = (su2double*)malloc(sizeof(su2double)*nBins);

        for(i=0;i<nBins;i++){
          diameters_multibin[i] = diameters_tmp[i];
          percentage_multibin[i] = percentage_tmp[i]; 
        }


        for(i=0;i<nBins;i++){
          diameters_multibin[i] = diameters_tmp[i];
          percentage_multibin[i] = percentage_tmp[i]; 
        }
      }
    }
    
    su2double LWC_inf = config_container[ZONE_0]->GetLiquidWaterContent();
    


    if(multibin){ //only steady multibin
 

      if(rank==MASTER_NODE){
        cout <<"+-----------------------------------------------------------------------+" << endl;
        cout << "|\n| Runnin particle tracking using following multibin data: \n|\n";
        cout <<"+-----------------------------------------------------------------------+" << endl;
        cout <<"| \t LWC (%) \t|     Particle Diameter (m)"<<endl;
        cout <<"+-----------------------------------------------------------------------+" << endl;
        for (i = 0;i<nBins;i++){
          cout << "| "<<i+1<<")\t " << percentage_multibin[i]<< "%   \t|     "<< diameters_multibin[i]<<"\n";
        
        }
        cout <<"+-----------------------------------------------------------------------+" << endl;
      }

      PTsolver->InitializeMultiBin(diameters_multibin, percentage_multibin, CFL_multibin, LWC_inf, nBins);

      //running simulation for each bin
      for (i = 0;i<nBins;i++){
        if(rank==MASTER_NODE){
          cout << "|\n|   Runnin Bin "<<i+1<<")\t" << percentage_multibin[i]<< "%   \t-    MVD = "<< diameters_multibin[i] << " (m)\n|\n";
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

            /*--- Sum this bin's collection efficiency ---*/
            PTsolver->AddBinCollectionEfficiency(geometry_fine);
            /*--- Compute BCs for splashing droplets (and corrects BETA) ---*/
            PTsolver->ComputeSplashingBCs(geometry_fine,config_container[ZONE_0],runtimeSplashing); 
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
        
        /*--- Compute BCs for splashing droplets (and corrects BETA) ---*/
        //cout<<"\n\n\t"<< rank <<" before ComputeSplashingBCs \n\n";
        PTsolver->ComputeSplashingBCs(geometry_fine,config_container[ZONE_0],runtimeSplashing);
          
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
    solver_container[ZONE_0][INST_0][MESH_0][PT_SOL]->SetInitialCondition(
    geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], TimeIter);
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
    iteration_container_PT[ZONE_0][INST_0]->Solve(output_container_PT[ZONE_0], integration_container, geometry_container,
                                               solver_container, numerics_container, config_container, surface_movement,
                                               grid_movement, FFDBox, ZONE_0, INST_0);  

    
  }

}

void CPTDriver::Postprocess() {

  if (runtimeFlow) {
    iteration_container[ZONE_0][INST_0]->Postprocess(
        output_container[ZONE_0], integration_container, geometry_container, solver_container, numerics_container,
        config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);
  } else {

    iteration_container_PT[ZONE_0][INST_0]->Postprocess(
        output_container_PT[ZONE_0], integration_container, geometry_container, solver_container, numerics_container,
        config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    
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
    iteration_container_PT[ZONE_0][INST_0]->Update(output_container_PT[ZONE_0], integration_container, geometry_container,
                                                solver_container, numerics_container, config_container,
                                                surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

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
    wrote_files = output_container_PT[ZONE_0]->SetResult_Files(geometry_container[ZONE_0][INST_0][MESH_0],
                                                               config_container[ZONE_0],
                                                               solver_container[ZONE_0][INST_0][MESH_0],
                                                               TimeIt, StopCalc);

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
    output = output_container_PT[ZONE_0];
    
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
    output = output_container_PT[ZONE_0];
    
  }
  return output->GetTimeConvergence();
}
