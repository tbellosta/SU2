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

/* GIUSEPPESIRIANNI */

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
/* GIUSEPPESIRIANNI */


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

    if(splashing){
      /*--- Set splashing as runtime ---*/
      setRuntimeSplashing();

      /*--- Compute BCs for splashing droplets ---*/
      ComputeSplashingBCs();
      

      cout << "\n before preprocess \n";
      Preprocess(TimeIter);

      /*--- Run a time-step iteration of the single-zone problem. ---*/
      cout << "\n before run \n";
      Run();

      /*--- Perform some postprocessing on the solution before the update ---*/
      cout << "\n before postprocess \n";
      Postprocess();

      /*--- Update the solution for dual time stepping strategy ---*/

      Update();

      /*--- Monitor the computations after each iteration. ---*/

      Monitor(TimeIter);

      /*--- Output the solution in files. ---*/
      Output(TimeIter);

    }
    


    /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    if (StopCalc) break;

    TimeIter++;

  }

}


//For now computes BCs for splashing droplets using Wright&Potacpzuk model
void CPTDriver::ComputeSplashingBCs() {

    /*--- Splashing Particle system ---*/
    unsigned short iDim, iVar;
    su2double mu_droplets = 18.03e-6;       //droplets fluid dynamic viscosity
    su2double rho_droplets = 1.0;           //droplets fluid density
    su2double sigma_droplets = 0.0756;      //droplet surface tension
    su2double diameter_splashing_droplets;  //splashing droplets diameter
    su2double diameter_droplets = 0;        //droplets diameter config->GetParticle_Size();
    su2double U_droplets;                   //droplets freestream droplets
    su2double *vinf = config_container[ZONE_0]->GetVelocity_FreeStream();
    
    U_droplets = GeometryToolbox::Norm(nDim, vinf);
    
    unsigned long iVertex, iPoint, Point_Normal;
    unsigned short iMarker, KindBC, val_marker;
    #define MAXNDIM 3
    su2double Area, ProjVelocity_i, V_boundary[4], V_domain[4], Normal[MAXNDIM] = {0.0}, UnitNormal[MAXNDIM] = {0.0},
                                                              V_wall[4], Relax;


    //Compute splashing droplets diameter (approx all splashing droplets of same diameter)
    //Compute Re_droplets
    su2double Re_droplets = 0;

    //Compute Ohnesorge number for splashing droplets
    su2double Oh_droplets = 0;
    
    //Compute K (mundo splashing parameter)
    su2double K = 0;
    su2double diameter_splashing_droplets_ALPHAAVG = 0;
    su2double tot_alpha = 0;



    //only finest mesh (?) single grid iteration (?)
    unsigned short FinestMesh = config_container[ZONE_0]->GetFinestMesh();    
    CGeometry* geometry_fine = geometry_container[ZONE_0][INST_0][FinestMesh];
    CSolver** solvers_fine = solver_container[ZONE_0][INST_0][FinestMesh];
    CSolver *splashingsolver = (solver_container[ZONE_0][INST_0][FinestMesh][SPLASHINGPT_SOL]);
    diameter_droplets = solvers_fine[PT_SOL]->GetDropletDiameter();
    
  
    //searching for euler wall marker (FOR NOW ONLY ONE EULER WALL ALLOWED (?))
    for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
      KindBC = config_container[ZONE_0]->GetMarker_All_KindBC(iMarker);
      if(KindBC == EULER_WALL){
        val_marker = iMarker;        
      }
    } 

    //Save splashing BC (ONLY ONE WALL ALLOWED AT THE MOMENT)
    
    //vector<vector<su2double>> splashingBCs(geometry_fine->nVertex[val_marker], vector<su2double>(3));//vector containing bcs for splashing droplets:   vec[ivertex][ivar]
    
    
    //Iterate on vertices of Wall BC ONLY 2D ONLY 2D ONLY 2D
    for (iVertex = 0; iVertex < geometry_fine->nVertex[val_marker]; iVertex++) {

      geometry_fine->vertex[val_marker][iVertex]->GetNormal(Normal);
      
      //for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Compute unit normal, to be used for unit tangential, projected velocity and velocity
            component gradients. ---*/
      Area = GeometryToolbox::Norm(nDim,Normal);

      for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = Normal[iDim] / Area;


      iPoint = geometry_fine->vertex[val_marker][iVertex]->GetNode();
      CVariable* nodes = solvers_fine[PT_SOL]->GetNodes();

      /*--- Get boundary solution at this boundary node of droplets---*/
      //check if not halo node
      if (geometry_fine->nodes->GetDomain(iPoint)) {
        unsigned short nVar = solvers_fine[PT_SOL]->GetnVar();
        for (iVar = 0; iVar < nVar; iVar++) V_domain[iVar] = nodes->GetPrimitive(iPoint,iVar);
      }

      U_droplets = sqrt(V_domain[1]*V_domain[1] + V_domain[2]*V_domain[2]);
      //Compute splashing droplets diameter (approx all splashing droplets of same diameter)
      //Compute Re_droplets
      su2double Re_droplets = (rho_droplets * U_droplets * diameter_droplets) / mu_droplets;

      //Compute Ohnesorge number for splashing droplets
      su2double Oh_droplets = mu_droplets / sqrt(rho_droplets * sigma_droplets * diameter_droplets);
      
      //Compute K (mundo splashing parameter)
      su2double K = Oh_droplets * pow(Re_droplets,1.25);
      diameter_splashing_droplets = diameter_droplets * 8.72 * exp(-0.0281 * K);
      //Compute for each vertex KW (LEWICE splashing parameter)
      su2double LWC = V_domain[0] * rho_droplets;
      su2double KW = pow(K,0.859) * pow((rho_droplets / LWC), 0.125);

      su2double theta;

      //Compute wall/droplets impact angle
      if(nDim==2){
        theta = 0.5*M_PI -  acos(abs(-V_domain[1]*UnitNormal[0]-V_domain[2]*UnitNormal[1]) / (sqrt(V_domain[1]*V_domain[1] + V_domain[2]*V_domain[2])));
      }
      else{
        theta = 0.5*M_PI - acos((-V_domain[1]*UnitNormal[0]-V_domain[2]*UnitNormal[1]-V_domain[3]*UnitNormal[2]) / (sqrt(V_domain[1]*V_domain[1] + V_domain[2]*V_domain[2] + V_domain[3]*V_domain[3])));
      }
      su2double theta_deg = (180 / M_PI * theta);
      
      su2double splashing_discriminator = 0;
      su2double theta_threshold = 0.001;          //threshold below which we will consider theta as basically 0 (no splashing can occur, floating point errors can occur)
      
      if(90 - abs(theta_deg) < theta_threshold){
        theta = 90;
        splashing_discriminator = -1;             //no splashing
      }else{
        splashing_discriminator = ( KW / pow(sin(abs(theta)), 1.25) ) - 200;
      }

      cout << "\ntheta = "<<theta_deg<<"deg\t Splash? = "<<(splashing_discriminator>0);
      cout << "\tdiameter = "<<diameter_splashing_droplets;

      su2double U_x_splashed = 0;
      su2double U_y_splashed = 0;
      su2double LWC_splashed = -1;
      if((diameter_splashing_droplets / diameter_droplets) < 1){
        if((diameter_splashing_droplets / diameter_droplets) > 0.05){
          //all ok
        }else{
          diameter_splashing_droplets = 0.05 * diameter_droplets;          
        }        
      }else{
          diameter_splashing_droplets = 1.0 * diameter_droplets;
      }

      //Verify if threshold is surpassed and therefore splashing occurs
      if(splashing_discriminator > 0){
        //splashing occurs
        su2double U_normal_domain = V_domain[1]*UnitNormal[0] + V_domain[2]*UnitNormal[1];
        su2double U_tangential_domain =  V_domain[1]*UnitNormal[1] - V_domain[2]*UnitNormal[0];
        LWC_splashed = LWC * 0.7 * (1 - sin(theta)) * (1 - exp(-0.0092026 * splashing_discriminator));
        su2double U_tangential_splashed = U_tangential_domain * (1.075 - 0.0025 * theta_deg);
        su2double U_normal_splashed = - U_normal_domain * (0.3 - 0.002 * theta_deg);
        U_x_splashed = U_normal_splashed * Normal[0] + U_tangential_splashed * Normal[1];
        U_y_splashed = U_normal_splashed * Normal[1] - U_tangential_splashed * Normal[0];


        tot_alpha += LWC_splashed / rho_droplets;
        diameter_splashing_droplets_ALPHAAVG += (LWC_splashed / rho_droplets)*diameter_splashing_droplets;
      }
      else{
        //no splashing occurs

      }

      //save in splashing solver the BCs
      splashingsolver->SetSplashingBCs(LWC_splashed / rho_droplets, U_x_splashed, U_y_splashed, iVertex);
  
      


      

    }

    splashingsolver->SetSplashingDiameter(diameter_splashing_droplets_ALPHAAVG / tot_alpha); 
      

    cout << "\n Diameter Splashing Droplets = "<<diameter_splashing_droplets<<"\n";
    cout << "\n Diameter Incoming Droplets = "<<diameter_droplets<<"\n";


    //All BCs are in a vector of vectors in the splashing solver container (?)


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

  cout << "\n before idef HAVE_MPI in CPTDriver::Preprocess() (splashing seems to have problems here)\n";

#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  cout << "\n after idef HAVE_MPI in CPTDriver::Preprocess()  \n";
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
