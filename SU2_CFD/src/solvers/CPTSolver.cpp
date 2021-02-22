/*!
 * \file CPTSolver.cpp
 * \brief Main subrotuines for solving the PT equations
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

#include "../../include/solvers/CPTSolver.hpp"

CPTSolver::CPTSolver(void) : CSolver() { }

CPTSolver::CPTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

  unsigned short iVar, iDim, nLineLets, iMarker;
  unsigned long iVertex;

  bool multizone = config->GetMultizone_Problem();

  int rank = MASTER_NODE;


  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();
  nMarker = config->GetnMarker_All();

  CurrentMesh = iMesh;

  /*--- Dimension of the problem  ---*/

  nVar = 3;
  if (nDim == 3) nVar++;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;
  nPrimVar = nVar;
  nPrimVarGrad = nVar;

  /*--- Define some auxiliary vector related with the residual ---*/

  Residual      = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_i    = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Residual_Max  = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Res_Conv      = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliar vector related with the solution ---*/

  Solution = new su2double[nVar];
  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];

  /*--- Define some auxiliary vectors related to the geometry ---*/

  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;

  /*--- Define some auxiliary vectors related to the primitive flow solution ---*/

  Primitive_i = new su2double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new su2double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) Primitive_j[iVar] = 0.0;

  /*--- Jacobians and vector structures for implicit computations ---*/

  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }

  if ((config->GetKind_ConvNumScheme_PT() == SPACE_CENTERED) && (MGLevel == MESH_0)) {
    iPoint_UndLapl = new su2double[nPointDomain];
    jPoint_UndLapl = new su2double[nPointDomain];
  }

  /*--- Initialization of the structure of the whole Jacobian ---*/

  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Eulerian Particle Tracking) MG level: " << iMesh << "." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

  if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  if (config->GetExtraOutput()) {
    if (nDim == 2) { nOutputVariables = 13; }          /** DON'T KNOW HOW MANY AND WHICH VARIABLES **/
    else if (nDim == 3) { nOutputVariables = 19; }
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
    OutputHeadingNames = new string[nOutputVariables];
  }

  /*--- Collection efficiency in all the markers ---*/

  CollectionEfficiency = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CollectionEfficiency[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CollectionEfficiency[iMarker][iVertex] = 0.0;
    }
  }

  if (multizone){
    /*--- Initialize the BGS residuals. ---*/
    Residual_BGS      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]  = 1.0;
    Residual_Max_BGS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
    }
  }

  /*--- Initialize the nodes vector. ---*/

  nodes = new CPTVariable(config->GetTemperature_FreeStreamND(), config->GetVelocity_FreeStreamND(), nPoint, nDim, nVar, config);

  SetBaseClassPointerToNodes();

  for (int iPoint = 0; iPoint < nPoint; ++iPoint) {
    nodes->SetLocalCFL(iPoint, config->GetCFL_PT());
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Add the solver name (max 8 characters) ---*/

  SolverName = "PT";

  V_inf = GeometryToolbox::Norm(nDim, config->GetVelocity_FreeStream());

  /*--- Check if we are executing a verification case. If so, the
 VerificationSolution object will be instantiated for a particular
 option from the available library of verification solutions. Note
 that this is done after SetNondim(), as problem-specific initial
 parameters are needed by the solution constructors. ---*/

  SetVerificationSolution(nDim, nVar, config);

}

CPTSolver::~CPTSolver(void) {

  unsigned short iMarker;

  if (CollectionEfficiency != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] CollectionEfficiency[iMarker];
    }
    delete [] CollectionEfficiency;
  }

  delete nodes;
}

void CPTSolver:: Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint;
  unsigned long InnerIter = config->GetInnerIter();

  bool center = (config->GetKind_ConvNumScheme_PT() == SPACE_CENTERED);
  bool muscl            = (config->GetMUSCL_PT());
  bool limiter          = (config->GetKind_SlopeLimit_PT() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool van_albada       = config->GetKind_SlopeLimit_PT() == VAN_ALBADA_EDGE;

//  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
//    /*--- Initialize the residual vector ---*/
//    LinSysRes.SetBlock_Zero(iPoint);
//  }
  LinSysRes.SetValZero();

  /*--- Initialize the Jacobian matrices ---*/
  Jacobian.SetValZero();

//  LimitSolution();

  /*--- Compute Primitives  ---*/
  SetPrimitiveVariables(geometry, config);


  if (config->GetReconstructionGradientRequired()) {
    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetPrimitive_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
  }
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetPrimitive_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetPrimitive_Gradient_LS(geometry, config);


  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    SetCentered_Dissipation_Sensor(geometry, config);
    SetUndivided_Laplacian(geometry, config);
  }

  if (limiter && (iMesh == MESH_0) && !Output && !van_albada) {

    SetPrimitive_Limiter(geometry, config);
  }


}

void CPTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  computeCollectionEfficiency(geometry,solver_container,config,iMesh);
}

void CPTSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/

  unsigned short iDim, iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;

  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;

  string restart_filename = config->GetFilename(config->GetSolution_FileName_PT(), "", val_iter);

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = 0;


  if (nDim == 2) skipVars = 2;
  if (nDim == 3) skipVars = 3;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local,Solution);
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    if (rank == MASTER_NODE) {
      cout << endl << "The solution file " << restart_filename.data() << " doesn't match with the mesh file!" << endl;
      cout << "It could be empty lines at the end of the file." << endl << endl;
    }
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We alo call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/

  solver[MESH_0][PT_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][PT_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  solver[MESH_0][PT_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_PT_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->nodes->GetVolume(iPoint);
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->nodes->GetnChildren_CV(iPoint); iChildren++) {
        Point_Fine = geometry[iMesh]->nodes->GetChildren_CV(iPoint, iChildren);
        Area_Children = geometry[iMesh-1]->nodes->GetVolume(Point_Fine);
        Solution_Fine = solver[iMesh-1][PT_SOL]->GetNodes()->GetSolution(Point_Fine);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][PT_SOL]->GetNodes()->SetSolution(iPoint,Solution);
    }
    solver[iMesh][PT_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][PT_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][PT_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_PT_SYS, false);
  }

  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars;
  delete [] Restart_Data;
  Restart_Vars = nullptr; Restart_Data = nullptr;

}


void CPTSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint, jPoint, iEdge;
  su2double *Diff;
  unsigned short iVar;
  bool boundary_i, boundary_j;

  Diff = new su2double[nVar];

  nodes->SetUnd_LaplZero();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    /*--- Solution differences ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = nodes->GetSolution(iPoint,iVar) - nodes->GetSolution(jPoint,iVar);

    boundary_i = geometry->nodes->GetPhysicalBoundary(iPoint);
    boundary_j = geometry->nodes->GetPhysicalBoundary(jPoint);

    /*--- Both points inside the domain, or both in the boundary ---*/

    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->nodes->GetDomain(iPoint)) nodes->SubtractUnd_Lapl(iPoint,Diff);
      if (geometry->nodes->GetDomain(jPoint)) nodes->AddUnd_Lapl(jPoint,Diff);
    }

    /*--- iPoint inside the domain, jPoint on the boundary ---*/

    if (!boundary_i && boundary_j)
      if (geometry->nodes->GetDomain(iPoint)) nodes->SubtractUnd_Lapl(iPoint,Diff);

    /*--- jPoint inside the domain, iPoint on the boundary ---*/

    if (boundary_i && !boundary_j)
      if (geometry->nodes->GetDomain(jPoint)) nodes->AddUnd_Lapl(jPoint,Diff);

  }

  /*--- MPI parallelization ---*/

  InitiateComms(geometry, config, UNDIVIDED_LAPLACIAN);
  CompleteComms(geometry, config, UNDIVIDED_LAPLACIAN);

  delete [] Diff;

}

void CPTSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container,  CNumerics **numerics_container,
                                    CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  CNumerics* numerics = numerics_container[CONV_TERM];

  su2double *V_i, *V_j, *W_i, *W_j;
  unsigned long iEdge, iPoint, jPoint;
  bool flow = ((config->GetKind_Solver() == INC_NAVIER_STOKES)
               || (config->GetKind_Solver() == INC_RANS)
               || (config->GetKind_Solver() == DISC_ADJ_INC_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_INC_RANS));


  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge ---*/
    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- set volumes to compute dissipation coefficients ---*/
    numerics->SetLaminarViscosity(geometry->nodes->GetVolume(iPoint),geometry->nodes->GetVolume(jPoint));


    /*--- Primitive variables w/o reconstruction ---*/
    V_i = nodes->GetPrimitive(iPoint);
    V_j = nodes->GetPrimitive(jPoint);

    W_i = nodes->GetSolution(iPoint);
    W_j = nodes->GetSolution(jPoint);

    numerics->SetUndivided_Laplacian(nodes->GetUndivided_Laplacian(iPoint), nodes->GetUndivided_Laplacian(jPoint));
    numerics->SetNeighbor(geometry->nodes->GetnNeighbor(iPoint), geometry->nodes->GetnNeighbor(jPoint));

    numerics->SetPrimitive(V_i, V_j);
    numerics->SetConservative(W_i, W_j);

    numerics->SetLambda(nodes->GetLambda(iPoint), nodes->GetLambda(jPoint));

    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

    LinSysRes.AddBlock(iPoint, Residual);
    LinSysRes.SubtractBlock(jPoint, Residual);

    /*--- Implicit part ---*/

    Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, Jacobian_i, Jacobian_j);
  }
}

void CPTSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                  CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[CONV_TERM];

  su2double *V_i, *V_j, *Sol_i, *Sol_i_Corrected, *Sol_j, *Sol_j_Corrected, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j,
      **Temp_i_Grad, **Temp_j_Grad, Project_Temp_i_Grad, Project_Temp_j_Grad;

  su2double MUSCLSol_i[4] = {0.0, 0.0, 0.0, 0.0};
  su2double MUSCLSol_j[4] = {0.0, 0.0, 0.0, 0.0};
  su2double Sol_zero[4] = {0.0, 0.0, 0.0, 0.0};


  unsigned short iDim, iVar;
  unsigned long iEdge, iPoint, jPoint;
  bool muscl = (config->GetMUSCL_PT());
  bool limiter = (config->GetKind_SlopeLimit_PT() != NO_LIMITER);
  bool van_albada = (config->GetKind_SlopeLimit_PT() == VAN_ALBADA_EDGE);

  su2double V_cent, V_upw_i, V_upw_j, *Limiter_i, *Limiter_j, DV_i, DV_j;

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      /*--- Points in edge ---*/
      iPoint = geometry->edges->GetNode(iEdge,0);
      jPoint = geometry->edges->GetNode(iEdge,1);
      numerics->SetNormal(geometry->edges->GetNormal(iEdge));

      V_i = nodes->GetPrimitive(iPoint);
      V_j = nodes->GetPrimitive(jPoint);

      /* Second order reconstruction */
      if (muscl) {

        for (iDim = 0; iDim < nDim; iDim++) {
          Vector_i[iDim] = 0.5*(geometry->nodes->GetCoord(jPoint, iDim) - geometry->nodes->GetCoord(iPoint, iDim));
        }

        Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
        Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

        if (limiter) {
          Limiter_i = nodes->GetLimiter_Primitive(iPoint);
          Limiter_j = nodes->GetLimiter_Primitive(jPoint);
        }

        /*Loop to correct the PT variables*/
        for (iVar = 0; iVar < nVar; iVar++) {

          /*Apply the Gradient to get the right solution value on the edge */
          Project_Grad_i = 0.0; Project_Grad_j = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
            Project_Grad_j += Vector_i[iDim]*Gradient_j[iVar][iDim];
          }

          /* Van Albada slope limiter */
          V_cent = V_j[iVar] - V_i[iVar];
//          V_upw_i = 2*Project_Grad_i - V_cent;
//          V_upw_j = 2*Project_Grad_j - V_cent;

          if (limiter) {
            /*--- beta values: 0 gradient extrapolation, -1 for MUSCL extrapolation, 0.5 3rd order extrapolation at face, 1/3 3rd order approximation to the derivative at the node ---*/
            const su2double beta = 0.0;
            if (van_albada) {

              //          Limiter_i = fmax(0.0, ( fabs(V_cent * V_upw_i) +  V_cent * V_upw_i + EPS) / (pow(V_cent,2) + pow(V_upw_i,2) + EPS));
              //          Limiter_j = fmax(0.0, ( fabs(V_cent * V_upw_j) +  V_cent * V_upw_j + EPS) / (pow(V_cent,2) + pow(V_upw_j,2) + EPS));

              Limiter_i[iVar] = V_cent * (2.0 * Project_Grad_i + V_cent) / (4 * pow(Project_Grad_i, 2) + pow(V_cent, 2) + EPS);
              Limiter_j[iVar] = V_cent * (2.0 * Project_Grad_j + V_cent) / (4 * pow(Project_Grad_j, 2) + pow(V_cent, 2) + EPS);
            }

//            MUSCLSol_i[iVar] = V_i[iVar] + 0.25 * Limiter_i[iVar] * ((1-beta*Limiter_i[iVar])*V_upw_i + (1+beta*Limiter_i[iVar])*V_cent);
//            MUSCLSol_j[iVar] = V_j[iVar] - 0.25 * Limiter_j[iVar] * ((1-beta*Limiter_j[iVar])*V_upw_j + (1+beta*Limiter_j[iVar])*V_cent);

            MUSCLSol_i[iVar] = V_i[iVar] + Limiter_i[iVar] * Project_Grad_i;
            MUSCLSol_j[iVar] = V_j[iVar] - Limiter_j[iVar] * Project_Grad_j;
          } else {
            MUSCLSol_i[iVar] = V_i[iVar] + Project_Grad_i;
            MUSCLSol_j[iVar] = V_j[iVar] - Project_Grad_j;
          }

        }

        /*--- Check for non-physical solutions after reconstruction. If found, use the
       cell-average value of the solution. This is a locally 1st order approximation,
       which is typically only active during the start-up of a calculation. ---*/

        bool neg_alpha_i = (MUSCLSol_i[0] < PT_EPS);
        bool neg_alpha_j = (MUSCLSol_j[0] < PT_EPS);

        nodes->SetNon_Physical(iPoint, neg_alpha_i);
        nodes->SetNon_Physical(jPoint, neg_alpha_j);


        /*--- Get updated state, in case the point recovered after the set. ---*/
//        neg_alpha_i = nodes->GetNon_Physical(iPoint);
//        neg_alpha_j = nodes->GetNon_Physical(jPoint);

        numerics->SetPrimitive(neg_alpha_i? V_i : MUSCLSol_i,  neg_alpha_j? V_j : MUSCLSol_j);

      }
      else {

        numerics->SetPrimitive(V_i, V_j);
      }

      numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      LinSysRes.AddBlock(iPoint, Residual);
      LinSysRes.SubtractBlock(jPoint, Residual);

      /*--- Implicit part ---*/

      Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, Jacobian_i, Jacobian_j);
    }

}

void CPTSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  return;

  CNumerics* numerics = numerics_container[VISC_TERM];

  su2double *V_i, *V_j, **V_i_Grad, **V_j_Grad;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;

  bool limiter = (config->GetKind_SlopeLimit_PT() != NO_LIMITER);
  bool van_albada = (config->GetKind_SlopeLimit_PT() == VAN_ALBADA_EDGE);


  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    /*--- Points coordinates, and normal vector ---*/

    numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                       geometry->nodes->GetCoord(jPoint));
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    V_i_Grad = nodes->GetGradient_Primitive(iPoint);
    V_j_Grad = nodes->GetGradient_Primitive(jPoint);
    numerics->SetPrimVarGradient(V_i_Grad, V_j_Grad);

    /*--- Primitive variables w/o reconstruction ---*/
    V_i = nodes->GetPrimitive(iPoint);
    V_j = nodes->GetPrimitive(jPoint);
    numerics->SetPrimitive(V_i, V_j);

    /*--- Compute numerical viscosity ---*/

    su2double dist = 0.0, Area = 0.0;
    for (iDim = 0; iDim < nDim; ++iDim) {
      dist += pow(geometry->nodes->GetCoord(jPoint,iDim) - geometry->nodes->GetCoord(iPoint,iDim),2);
      Area += pow(geometry->edges->GetNormal(iEdge)[iDim],2);
    }
    Area = sqrt(Area);
    dist = 0.5*sqrt(dist);

    su2double impulse_two, impulse_single, impulse, Visc_i, Visc_j;
    su2double u_star, Proj_vel_i, Proj_vel_j, UnitNormal[3];

    Proj_vel_i = Proj_vel_j = 0.0;

    for (iDim = 0; iDim < nDim; ++iDim) {
      UnitNormal[iDim] = geometry->edges->GetNormal(iEdge)[iDim]/Area;
    }

    for (iDim = 0; iDim < nDim; ++iDim) {
      Proj_vel_i += V_i[iDim+1]*UnitNormal[iDim];
      Proj_vel_j += V_j[iDim+1]*UnitNormal[iDim];
    }

    u_star = (V_i[0]*Proj_vel_i + V_j[0]*Proj_vel_j) / (V_i[0] + V_j[0] + EPS);

    impulse_two = V_i[0]*Proj_vel_i - V_j[0]*Proj_vel_j;
    impulse_single = (V_i[0] - V_j[0])*u_star;

    impulse = fmax(0,impulse_two-impulse_single);

    Visc_i = impulse * dist / ((V_i[0]+V_j[0]) - V_j[0] + EPS);
    Visc_j = impulse * dist / ((V_i[0]+V_j[0]) - V_i[0] + EPS);

//    numerics->SetLaminarViscosity(1,1);
    numerics->SetLaminarViscosity(0.002,0.002);

    /*--- Compute residual, and Jacobians ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

    /*--- Add and subtract residual, and update Jacobians ---*/

//    LinSysRes.SubtractBlock(iPoint, Residual);
//    LinSysRes.AddBlock(jPoint, Residual);
//
//    Jacobian.UpdateBlocksSub(iEdge, iPoint, jPoint, Jacobian_i, Jacobian_j);


  }

}


void CPTSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                               unsigned short iMesh, unsigned long Iteration) {

  const bool implicit      = (config->GetKind_TimeIntScheme_PT() == EULER_IMPLICIT);
  const bool time_stepping = (config->GetTime_Marching() == TIME_STEPPING);
  const bool dual_time     = (config->GetTime_Marching() == DT_STEPPING_1ST) ||
                             (config->GetTime_Marching() == DT_STEPPING_2ND);

  su2double Global_Delta_UnstTimeND = 0.0;

  /*--- Init thread-shared variables to compute min/max values.
   *    Critical sections are used for this instead of reduction
   *    clauses for compatibility with OpenMP 2.0 (Windows...). ---*/


  Min_Delta_Time = 1e30;
  Max_Delta_Time = 0.0;
  Global_Delta_UnstTimeND = 1e30;


  const su2double *Normal = nullptr;
  su2double Area, Vol, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time, Local_Delta_Time_Visc;
  su2double Mean_LaminarVisc, Mean_EddyVisc, Mean_Density, Lambda_1, Lambda_2, Unitnormal[3];
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

  /*--- Loop domain points. ---*/

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    /*--- Set maximum eigenvalues to zero. ---*/

    nodes->SetMax_Lambda_Inv(iPoint,0.0);

    /*--- Loop over the neighbors of point i. ---*/

    for (unsigned short iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); ++iNeigh) {
      jPoint = geometry->nodes->GetPoint(iPoint,iNeigh);

      iEdge = geometry->nodes->GetEdge(iPoint,iNeigh);
      Normal = geometry->edges->GetNormal(iEdge);
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += pow(Normal[iDim],2); Area = sqrt(Area);

      /*--- Mean Values ---*/

      Mean_ProjVel = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));

      /*--- Adjustment for grid movement ---*/

      if (dynamic_grid) {
        const su2double *GridVel_i = geometry->nodes->GetGridVel(iPoint);
        const su2double *GridVel_j = geometry->nodes->GetGridVel(jPoint);

        for (iDim = 0; iDim < nDim; iDim++)
          Mean_ProjVel -= 0.5 * (GridVel_i[iDim] + GridVel_j[iDim]) * Normal[iDim];
      }

      /*--- Inviscid contribution ---*/

      Lambda = fabs(Mean_ProjVel);
      nodes->AddMax_Lambda_Inv(iPoint,Lambda);
    }

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        /*--- Point identification, Normal vector and area ---*/

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (!geometry->nodes->GetDomain(iPoint)) continue;
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

        /*--- Mean Values ---*/

        Mean_ProjVel = nodes->GetProjVel(iPoint,Normal);

        /*--- Adjustment for grid movement ---*/

        if (dynamic_grid) {
          const su2double *GridVel = geometry->nodes->GetGridVel(iPoint);

          for (iDim = 0; iDim < nDim; iDim++)
            Mean_ProjVel -= GridVel[iDim]*Normal[iDim];
        }

        /*--- Inviscid contribution ---*/

        Lambda = fabs(Mean_ProjVel);
        nodes->AddMax_Lambda_Inv(iPoint,Lambda);

      }
    }
  }

  /*--- Each element uses their own speed, steady state simulation. ---*/
  Min_Delta_Time = 1.0e20;
  Max_Delta_Time = 0.0;
  su2double Global_Delta_Time = 0.0;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->nodes->GetVolume(iPoint);

    if (Vol != 0.0) {
      Local_Delta_Time = nodes->GetLocalCFL(iPoint)*Vol / nodes->GetMax_Lambda_Inv(iPoint);

      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);

      nodes->SetDelta_Time(iPoint, min(Local_Delta_Time, config->GetMax_DeltaTime()));
    }
    else {
      nodes->SetDelta_Time(iPoint,0.0);
    }
  }

  /*--- Compute the min/max dt (in parallel, now over mpi ranks). ---*/

  if (config->GetComm_Level() == COMM_FULL) {
    su2double rbuf_time;
    SU2_MPI::Allreduce(&Min_Delta_Time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;

    SU2_MPI::Allreduce(&Max_Delta_Time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
  }

  config->SetDelta_UnstTimeND(config->GetDelta_UnstTime());

//  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
//    nodes->SetDelta_Time(iPoint,Min_Delta_Time);
//  }

  /*--- For exact time solution use the minimum delta time of the whole mesh. ---*/
  if (time_stepping) {

    /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step,
     *    otherwise it computes the time step based on the unsteady CFL. ---*/

    if (config->GetUnst_CFL() == 0.0) {
      Global_Delta_Time = config->GetDelta_UnstTime();
    }
    else {
      Global_Delta_Time = Min_Delta_Time;
    }
    Max_Delta_Time = Global_Delta_Time;

    config->SetDelta_UnstTimeND(Global_Delta_Time);

    /*--- Sets the regular CFL equal to the unsteady CFL. ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      nodes->SetLocalCFL(iPoint, config->GetUnst_CFL());
      nodes->SetDelta_Time(iPoint, Global_Delta_Time);
    }

  }

  /*--- Recompute the unsteady time step for the dual time strategy if the unsteady CFL is diferent from 0. ---*/

  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {

    /*--- Thread-local variable for reduction. ---*/
    su2double glbDtND = 1e30;
    Global_Delta_UnstTimeND = 1e30;

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Global_Delta_UnstTimeND = min(Global_Delta_UnstTimeND, config->GetUnst_CFL()*Global_Delta_Time / nodes->GetLocalCFL(iPoint));
    }

    SU2_MPI::Allreduce(&Global_Delta_UnstTimeND, &glbDtND, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = glbDtND;

    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }

  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/

  if (dual_time && !implicit) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), nodes->GetDelta_Time(iPoint));
      nodes->SetDelta_Time(iPoint, Local_Delta_Time);
    }
  }

}

void CPTSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res, alpha, deltaAlpha, factor;
  unsigned short iVar;
  unsigned long iPoint;

  bool adjoint = config->GetContinuous_Adjoint();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);;
    Delta = nodes->GetDelta_Time(iPoint) / Vol;

    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    local_Residual = LinSysRes.GetBlock(iPoint);

    /*--- Limit solution update to enforce positive volume fraction ---*/

    alpha = nodes->GetSolution(iPoint,0);
    Res = local_Residual[0] + local_Res_TruncError[0];
    Res = -Res*Delta;
    deltaAlpha = max(Res, PT_EPS-alpha);
    factor = (Res) ? deltaAlpha / Res : 1.0;

//    factor = 1.0;

//    if (alpha+Res <= PT_EPS) {throw std::runtime_error("error");}

    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = local_Residual[iVar] + local_Res_TruncError[iVar];
        Res = -Res*Delta;
        nodes->AddSolution(iPoint,iVar, Res*factor);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
      }
    }

  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

  /*--- For verification cases, compute the global error metrics. ---*/

  ComputeVerificationError(geometry, config);

}


void CPTSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, Vol, *local_Res_TruncError, **jacobian, tau;


  /*--- Set maximum residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- allocate jacobian ---*/
  jacobian = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    jacobian[iVar] = new su2double[nVar];

  for (iVar = 0; iVar < nVar; ++iVar)
    for (unsigned short jVar = 0; jVar < nVar; ++jVar)
      jacobian[iVar][jVar] = 0.0;

  /*--- Build implicit system ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the residual ---*/

    local_Res_TruncError = nodes->GetResTruncError(iPoint);

    /*--- Read the volume ---*/

    Vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);

    /*--- Modify matrix diagonal to assure diagonal dominance ---*/

    if (nodes->GetDelta_Time(iPoint) != 0.0) {

      Delta = Vol / nodes->GetDelta_Time(iPoint);

//      tau = computeRelaxationTime(solver_container,iPoint);
//      if (tau > nodes->GetDelta_Time(iPoint)) {
//        Jacobian.AddVal2Diag(iPoint, Delta);
//      } else {
//
//        jacobian[0][0] = Delta;
//        for (iVar = 1; iVar < nVar; ++iVar) {
//          jacobian[iVar][iVar] = Vol / max(1.0e-8, tau);
//        }
//        Jacobian.AddBlock2Diag(iPoint, jacobian);
//
//      }

      Jacobian.AddVal2Diag(iPoint, Delta);

    } else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index]);
//      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
    }
  }

  /*--- Initialize residual and solution at the ghost points ---*/

  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }

  /*--- Solve or smooth the linear system ---*/

  auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  SetIterLinSolver(iter);
  SetResLinSolver(System.GetResidual());

  /*--- Limit solution update to enforce positive volume fraction ---*/

//  LimitSolution();

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    for (iVar = 0; iVar < nVar; iVar++)
      nodes->AddSolution(iPoint,iVar, LinSysSol[iPoint*nVar+iVar]);

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

  /*--- For verification cases, compute the global error metrics. ---*/

  ComputeVerificationError(geometry, config);

  for (iVar = 0; iVar < nVar; ++iVar) delete [] jacobian[iVar];
  delete [] jacobian;

}

void CPTSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar;
  su2double Area_Children, Area_Parent, *Solution_Fine, *Solution;

  bool restart   = (config->GetRestart() && config->GetRestart_PT());
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));


  if (!restart && (TimeIter == 0)) {

    su2double* Coord;
    su2double LWC = config->GetLiquidWaterContent();

    su2double *flowPrimitive;


    for (iMesh = 0; iMesh <= config->GetnMGLevels(); ++iMesh) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); ++iPoint) {

        Coord = geometry[iMesh]->nodes->GetCoord(iPoint);
        su2double *solDOF = nodes->GetSolution(iPoint);

        if (VerificationSolution) {

          /* Set the solution in this DOF to the initial condition provided by
             the verification solution class. This can be the exact solution,
             but this is not necessary. */
          VerificationSolution->GetInitialCondition(Coord, solDOF);
          VerificationSolution->GetPrimitive(Coord, 0.0, nodes->GetPrimitive(iPoint));

        } else {
          LWC = 1.0;
          solDOF[0] = LWC;

          for (int iDim = 0; iDim < nDim; ++iDim) solDOF[iDim+1] = 1 * config->GetVelocity_FreeStream()[iDim];

//        flowPrimitive = solver_container[iMesh][FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
//        for (int iDim = 0; iDim < nDim; ++iDim) solDOF[iDim+1] = flowPrimitive[iDim+1] * LWC;
        }

      }
//      solver_container[iMesh][PT_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
//      solver_container[iMesh][PT_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    }
  }

  /*--- The value of the solution for the first iteration of the dual time ---*/

  if (dual_time && (TimeIter == 0 || (restart && (long)TimeIter == (long)config->GetRestart_Iter()))) {

    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/

    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      solver_container[iMesh][PT_SOL]->GetNodes()->Set_Solution_time_n();
      solver_container[iMesh][PT_SOL]->GetNodes()->Set_Solution_time_n1();
    }

    if ((restart && (long)TimeIter == (long)config->GetRestart_Iter()) &&
        (config->GetTime_Marching() == DT_STEPPING_2ND)) {

      /*--- Load an additional restart file for a 2nd-order restart ---*/

      solver_container[MESH_0][PT_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetRestart_Iter()-1), true);

      /*--- Push back this new solution to time level N. ---*/

      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        solver_container[iMesh][PT_SOL]->GetNodes()->Set_Solution_time_n();
      }
    }
  }
}

void CPTSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                       unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {

  /*--- Local variables ---*/

  unsigned short iVar, jVar;
  unsigned long iPoint;

  su2double *U_time_n, *U_time_nP1, *U_time_nM1;
  su2double Volume_nP1, TimeStep;

  bool implicit       = (config->GetKind_TimeIntScheme_PT() == EULER_IMPLICIT);

  /*--- Store the physical time step ---*/

  TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the dual time-stepping source term for static meshes ---*/

  if (!dynamic_grid) {

    /*--- Loop over all nodes (excluding halos) ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/

      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n   = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);

      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/

      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetTime_Marching() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetTime_Marching() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                             +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }

      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/

      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetTime_Marching() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetTime_Marching() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }

        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
      }
    }
  }
}


void CPTSolver::BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                             CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {

  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  su2double *V_boundary, *V_domain, V_Infty[4], Vn;

  bool implicit             = (config->GetKind_TimeIntScheme_PT() == EULER_IMPLICIT);

  su2double Normal[3], LWC = 1.0;

  V_Infty[0] = LWC;
  for (iDim = 0; iDim < nDim; ++iDim) {
    V_Infty[iDim+1] = config->GetVelocity_FreeStream()[iDim];
  }

  su2double time = config->GetPhysicalTime();


  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->nodes->GetDomain(iPoint)) {

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      conv_numerics->SetNormal(Normal);

      Vn = 0;
      for (iDim = 0; iDim < nDim; ++iDim) {
        Vn -= nodes->GetPrimitive(iPoint,iDim+1)*Normal[iDim];
      }

      if (VerificationSolution) {
        VerificationSolution->GetPrimitive(geometry->nodes->GetCoord(iPoint), time, V_Infty);
      }

      /*--- Retrieve solution at this boundary node ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      if (Vn >= 0) {
        V_boundary = V_Infty;
      } else V_boundary = nodes->GetPrimitive(iPoint);

//      su2double VV[3];
//      VV[0] = V_Infty[0];
//      VV[1] = uEx(geometry->nodes->GetCoord(iPoint));
//      VV[2] = vEx(geometry->nodes->GetCoord(iPoint));
//
//      nodes->SetSolution_Old(iPoint, V_Infty);
//      nodes->SetSolution(iPoint, V_Infty);
//      nodes->SetPrimitive(iPoint, VV);
//      nodes->SetRes_TruncErrorZero(iPoint);
//      LinSysRes.SetBlock_Zero(iPoint);
//
//      /*--- Adjust rows of the Jacobian (includes 1 in the diagonal) ---*/
//
//      if (implicit) {
//        for (iVar = 0; iVar < nVar; iVar++) {
//          total_index = iPoint * nVar + iVar;
//          Jacobian.DeleteValsRowi(total_index);
//        }
//      }
//
//      return;


      conv_numerics->SetPrimitive(nodes->GetPrimitive(iPoint), V_boundary);

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

      /*---- viscous part ----*/

      /*--- Set the normal vector and the coordinates ---*/

//      visc_numerics->SetNormal(Normal);
//      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(Point_Normal));
//
//      /*--- Primitive variables, and gradient ---*/
//
//      visc_numerics->SetPrimitive(nodes->GetPrimitive(iPoint), V_boundary);
//      visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
//                                        nodes->GetGradient_Primitive(iPoint));
//
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Update residual value ---*/
//      LinSysRes.SubtractBlock(iPoint, Residual);
//
//      /*--- Jacobian contribution for implicit integration. ---*/
//      if (implicit)
//        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);

    }
  }

}
void CPTSolver::Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics_container,
                                CConfig* config, unsigned short iMesh) {

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM];

  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  su2double tau;

  /*--- Loop over all points. ---*/
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    numerics->SetConservative(nodes->GetPrimitive(iPoint), nullptr);

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    tau = computeRelaxationTime(solver_container, iPoint);
    numerics->SetParticleTau(tau);

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Compute the source term ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, Residual);

    Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);

  }

  /*--- Check if a verification solution is to be computed. ---*/

  if ( VerificationSolution ) {
    if ( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Get the physical time. ---*/
      su2double time = 0.0;
      if (config->GetTime_Marching()) time = config->GetPhysicalTime();

      /*--- Loop over points ---*/
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

        /*--- Get control volume size. ---*/
        su2double Volume = geometry->nodes->GetVolume(iPoint);

        /*--- Get the current point coordinates. ---*/
        const su2double *coor = geometry->nodes->GetCoord(iPoint);

        /*--- Get the MMS source term. ---*/
        vector<su2double> sourceMan(nVar,0.0);
        VerificationSolution->GetMMSSourceTerm(coor, time, sourceMan.data());

        /*--- Compute the residual for this control volume and subtract. ---*/
        for (unsigned long iVar = 0; iVar < nVar; iVar++) {
          LinSysRes(iPoint,iVar) -= sourceMan[iVar]*Volume;
        }
      }
    }
  }

}


void CPTSolver::BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                            CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint, Point_Normal;

#define MAXNDIM 3

  bool implicit = (config->GetKind_TimeIntScheme_PT() == EULER_IMPLICIT);
  bool viscous = config->GetViscous();
  bool preprocessed = false;

  /*--- Allocation of variables necessary for convective fluxes. ---*/
  su2double Area, ProjVelocity_i, V_reflected[4], *V_domain, Normal[MAXNDIM] = {0.0}, UnitNormal[MAXNDIM] = {0.0};

  /*--- Loop over all the vertices on this boundary marker. ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    if (!preprocessed || !geometry->bound_is_straight[val_marker]) {
      /*----------------------------------------------------------------------------------------------*/
      /*--- Preprocessing:                                                                         ---*/
      /*--- Compute the unit normal and (in case of viscous flow) a corresponding unit tangential  ---*/
      /*--- to that normal. On a straight(2D)/plane(3D) boundary these two vectors are constant.   ---*/
      /*--- This circumstance is checked in gemoetry->ComputeSurf_Straightness(...) and stored     ---*/
      /*--- such that the recomputation does not occur for each node. On true symmetry planes, the ---*/
      /*--- normal is constant but this routines is used for Symmetry, Euler-Wall in inviscid flow ---*/
      /*--- and Euler Wall in viscous flow as well. In the latter curvy boundaries are likely to   ---*/
      /*--- happen. In doubt, the conditional above which checks straightness can be thrown out    ---*/
      /*--- such that the recomputation is done for each node (which comes with a tiny performance ---*/
      /*--- penalty).                                                                              ---*/
      /*----------------------------------------------------------------------------------------------*/

      preprocessed = true;

      /*--- Normal vector for a random vertex (zero) on this marker (negate for outward convention). ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Compute unit normal, to be used for unit tangential, projected velocity and velocity
            component gradients. ---*/
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim] * Normal[iDim];
      Area = sqrt(Area);

      for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim] / Area;

    }      // if bound_is_straight

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {
      /*-------------------------------------------------------------------------------*/
      /*--- Step 1: For the convective fluxes, create a reflected state of the      ---*/
      /*---         Primitive variables by copying all interior values to the       ---*/
      /*---         reflected. Only the velocity is mirrored along the symmetry     ---*/
      /*---         axis. Based on the Upwind_Residual routine.                     ---*/
      /*-------------------------------------------------------------------------------*/

//      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Grid movement ---*/
      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      /*--- Normal vector for this vertex (negate for outward convention). ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Get current solution at this boundary node ---*/
      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Set the reflected state based on the boundary node. Scalars are copied and
            the velocity is mirrored along the symmetry boundary, i.e. the velocity in
            normal direction is substracted twice. ---*/
      for (iVar = 0; iVar < nVar; iVar++) V_reflected[iVar] = V_domain[iVar];

      /*--- Compute velocity in normal direction (ProjVelcity_i=(v*n)) und substract twice from
            velocity in normal direction: v_r = v - 2 (v*n)n ---*/
      ProjVelocity_i = nodes->GetProjVel(iPoint,UnitNormal);

      /*--- Adjustment to v.n due to grid movement. ---*/
      if (dynamic_grid) {
        ProjVelocity_i -= GeometryToolbox::DotProduct(nDim, geometry->nodes->GetGridVel(iPoint), UnitNormal);
      }

      if (ProjVelocity_i >= 0) {
        for (iDim = 0; iDim < nDim; iDim++)
          V_reflected[iDim + 1] = nodes->GetPrimitive(iPoint, iDim + 1) - 2.0 * ProjVelocity_i * UnitNormal[iDim];
      }

      /*--- Set Primitive and Secondary for numerics class. ---*/
      conv_numerics->SetPrimitive(V_domain, V_reflected);

      /*--- Compute the residual using an upwind scheme. ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration. ---*/
      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);



      /*---- viscous part ----*/

      /*--- Set the normal vector and the coordinates ---*/

//      visc_numerics->SetNormal(Normal);
//      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(Point_Normal));
//
//      /*--- Primitive variables, and gradient ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_reflected);
//      visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
//                                        nodes->GetGradient_Primitive(iPoint));
//
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Update residual value ---*/
//      LinSysRes.SubtractBlock(iPoint, Residual);
//
//      /*--- Jacobian contribution for implicit integration. ---*/
//      if (implicit)
//        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);


    }    // if GetDomain
  }      // for iVertex

}


void CPTSolver::computeCollectionEfficiency(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CConfig *config,
                                            unsigned short iMesh) {

  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint;

  su2double Normal[3], Area, NDfactor;

  NDfactor = V_inf;

  for (iMarker = 0; iMarker < nMarker; ++iMarker) {
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; ++iVertex) {

      geometry->vertex[iMarker][iVertex]->GetNormal(Normal);

      /*--- Compute unit normal ---*/
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim] * Normal[iDim];
      Area = sqrt(Area);

      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim] / Area;

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

      CollectionEfficiency[iMarker][iVertex] = 0.0;

      for (iDim = 0; iDim < nDim; ++iDim)
//        CollectionEfficiency[iMarker][iVertex] += nodes->GetPrimitive(iPoint, iDim + 1) * Normal[iDim];
        CollectionEfficiency[iMarker][iVertex] += nodes->GetSolution(iPoint, iDim + 1) * Normal[iDim];

      CollectionEfficiency[iMarker][iVertex] /= NDfactor;

    }
  }
}

void CPTSolver::BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                         CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry,solver_container,conv_numerics,visc_numerics,config,val_marker);

}

void CPTSolver::LimitSolution(void) {

  unsigned long iPoint;
  unsigned short iDim, iVar;

  su2double alpha, deltaAlpha, factor;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    alpha = nodes->GetSolution(iPoint,0);
    deltaAlpha = max(LinSysSol[iPoint*nVar], PT_EPS-alpha);
//    deltaAlpha = min(max(LinSysSol[iPoint*nVar], PT_EPS-alpha), 2*alpha);
    factor = (LinSysSol[iPoint*nVar]) ? deltaAlpha / LinSysSol[iPoint*nVar] : 1.0;

//    if (alpha+LinSysSol[iPoint*nVar] <= PT_EPS) {throw std::runtime_error("error");}

    LinSysSol[iPoint*nVar] = deltaAlpha;

    for (iVar = 1; iVar < nVar; iVar++) LinSysSol[iPoint*nVar+iVar] = LinSysSol[iPoint*nVar+iVar]*factor;
  }

}
void CPTSolver::SetPrimitiveVariables(CGeometry *geometry, CConfig *config) {

  su2double alpha;
  su2double Sol_zero[4] = {0.0, 0.0, 0.0, 0.0};

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    alpha = nodes->GetSolution(iPoint,0);
    nodes->SetPrimitive(iPoint, 0, alpha);

    if (alpha != 0) {
      for (int iVar = 1; iVar < nPrimVar; ++iVar) {
        nodes->SetPrimitive(iPoint, iVar, nodes->GetSolution(iPoint, iVar) / alpha);
      }
    } else {
      for (int iVar = 0; iVar < nPrimVar; ++iVar) {
//        nodes->SetPrimitive(iPoint, iVar, nodes->GetSolution(iPoint, iVar));
        nodes->SetPrimitive(iPoint, iVar, Sol_zero[iVar]);
//        nodes->SetSolution(iPoint, iVar, Sol_zero[iVar]);
//        nodes->SetSolution_Old(iPoint, iVar, Sol_zero[iVar]);
      }
    }
  }

}

void CPTSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, const CConfig *config) {


  /*--- We can access memory more efficiently if there are no periodic boundaries. ---*/

  const bool isPeriodic = (config->GetnMarker_Periodic() > 0);

  /*--- Loop domain points. ---*/

  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    const bool boundary_i = geometry->nodes->GetPhysicalBoundary(iPoint);
    const su2double Density_i = nodes->GetPrimitive(iPoint,0);

    /*--- Initialize. ---*/
    iPoint_UndLapl[iPoint] = 0.0;
    jPoint_UndLapl[iPoint] = 0.0;

    /*--- Loop over the neighbors of point i. ---*/
    for (auto jPoint : geometry->nodes->GetPoints(iPoint))
    {
      bool boundary_j = geometry->nodes->GetPhysicalBoundary(jPoint);

      /*--- If iPoint is boundary it only takes contributions from other boundary points. ---*/
      if (boundary_i && !boundary_j) continue;

      su2double Density_j = nodes->GetPrimitive(jPoint,0);

      /*--- Dissipation sensor, add pressure difference and pressure sum. ---*/
      iPoint_UndLapl[iPoint] += Density_j - Density_i;
      jPoint_UndLapl[iPoint] += Density_j + Density_i;
    }

    if (!isPeriodic) {
      nodes->SetSensor(iPoint, fabs(iPoint_UndLapl[iPoint]) / (4*Density_i));
    }
  }

  if (isPeriodic) {
    /*--- Correct the sensor values across any periodic boundaries. ---*/

    for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
      InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_SENSOR);
      CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_SENSOR);
    }

    /*--- Set final pressure switch for each point ---*/

    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
      nodes->SetSensor(iPoint, fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);
  }

  /*--- MPI parallelization ---*/

  InitiateComms(geometry, config, SENSOR);
  CompleteComms(geometry, config, SENSOR);

}

void CPTSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {

  /*--- Loop domain points. ---*/

  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    /*--- Set eigenvalues to zero. ---*/
    nodes->SetLambda(iPoint,0.0);

    /*--- Loop over the neighbors of point i. ---*/
    for (unsigned short iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); ++iNeigh)
    {
      auto jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);

      auto iEdge = geometry->nodes->GetEdge(iPoint, iNeigh);
      auto Normal = geometry->edges->GetNormal(iEdge);
      su2double Area = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; iDim++) Area += pow(Normal[iDim],2);
      Area = sqrt(Area);

      /*--- Mean Values ---*/

      su2double Mean_ProjVel = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));

      /*--- Adjustment for grid movement ---*/

      if (dynamic_grid) {
        const su2double *GridVel_i = geometry->nodes->GetGridVel(iPoint);
        const su2double *GridVel_j = geometry->nodes->GetGridVel(jPoint);

        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          Mean_ProjVel -= 0.5 * (GridVel_i[iDim] + GridVel_j[iDim]) * Normal[iDim];
      }

      /*--- Inviscid contribution ---*/

      su2double Lambda = fabs(Mean_ProjVel);
      nodes->AddLambda(iPoint, Lambda);
    }

  }

  /*--- Loop boundary edges ---*/

  for (unsigned short iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

      for (unsigned long iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        /*--- Point identification, Normal vector and area ---*/

        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        su2double Area = 0.0;
        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          Area += pow(Normal[iDim],2);
        Area = sqrt(Area);

        /*--- Mean Values ---*/

        su2double Mean_ProjVel = nodes->GetProjVel(iPoint,Normal);

        /*--- Adjustment for grid movement ---*/

        if (dynamic_grid) {
          auto GridVel = geometry->nodes->GetGridVel(iPoint);
          for (unsigned short iDim = 0; iDim < nDim; iDim++)
            Mean_ProjVel -= GridVel[iDim]*Normal[iDim];
        }

        /*--- Inviscid contribution ---*/

        su2double Lambda = fabs(Mean_ProjVel);
        if (geometry->nodes->GetDomain(iPoint)) {
          nodes->AddLambda(iPoint,Lambda);
        }
      }
    }
  }

  /*--- Correct the eigenvalue values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_MAX_EIG);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_MAX_EIG);
  }

  /*--- MPI parallelization ---*/

  InitiateComms(geometry, config, MAX_EIGENVALUE);
  CompleteComms(geometry, config, MAX_EIGENVALUE);

}

su2double CPTSolver::computeRelaxationTime(CSolver** solver_container, unsigned long iPoint) {

  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  const su2double rho = 1000;
  const su2double sigma = 0.0756;
//  const su2double d = config->GetParticle_Size();
  const su2double d = 2e-5;
  su2double mu = 18.03e-6;
  su2double *FlowPrim = flowNodes->GetPrimitive(iPoint);
  su2double *Prim = nodes->GetPrimitive(iPoint);

  su2double *uFlow = &FlowPrim[1], rhoFlow = FlowPrim[nDim+2];

//  rhoFlow = 1.2;
//  su2double Uf[2] = {4.0, 1.0};
//  uFlow = Uf;

  su2double V_mag = 0.0;
  for (int iDim = 0; iDim < nDim; ++iDim) {
    V_mag += pow(uFlow[iDim]-Prim[iDim+1],2);
  }
  V_mag = sqrt(V_mag);

  su2double Re = rhoFlow*d*V_mag / mu;
//  Re = min(1.0e6,Re);
  su2double Cd;

  const su2double We = d * pow(V_mag,2) * rhoFlow / sigma;

  /*---------
  *  The function use the Morris approximation up to Re = 1e6 and then
  *  it is linked to the Shankar Subramanian approximation for Re > 1e6.
  *
  * Extended drag model accounting for the increase of drag
  * due to droplet deformation in SLD conditions
  * the model is based on droplet eccentricity and
  * on the drag coeficient of the sphere and the disk
  -----------*/


  /*---- Drag coefficient of the sphere (Morris+Shankar Subramanian) ----*/

  su2double Cd_sphere;

  const su2double delta = 0.0315;   /*---- correction factor to link the two curves ----*/

  if ( Re == 0 )
    Cd_sphere = 1;     /*---- any value is ok because the resulting drag is zero ----*/
  else if ( Re < 1.0e6)
    Cd_sphere =  24.0 / Re +
                 2.6 * (Re / 5.0) / (1 + pow((Re / 5.0), 1.52) ) +
                 0.411 * pow((Re / 263000.0), -7.94) / (1 + pow((Re / 263000.0), -8.0) ) +
                 pow(Re, 0.8) / 461000.0;
  else if ( Re > 1.0e6 )
    Cd_sphere = 0.19 - 8.0e4 / Re + delta;


  /*---- Drag coefficient of the disk (Clift) ----*/
  su2double Cd_disk;
  su2double x;

  if ( Re == 0 )
    Cd_disk = 1;     /*---- any value is ok because the resulting drag is zero  ----*/
  else if ( Re <= 0.01)
    Cd_disk =  64.0 / (M_PI*Re) *(1 + (Re/(2*M_PI)));
  else if ( Re > 0.01 && Re <= 1.5 )
  {
    x = -0.883 + 0.906*log10(Re)-0.025*pow(log10(Re),2);
    Cd_disk = 64.0 / (M_PI*Re) *(1 + pow(10,x));
  }
  else if ( Re > 1.5 && Re <= 133 )
    Cd_disk = 64.0 / (M_PI*Re) *(1 + 0.138*pow(Re,0.792));
  else if ( Re > 133)
    Cd_disk = 1.17;


  /*---- Eccentricity function ----*/
  su2double eccentricity;
  eccentricity = 1.0 - pow(1.0 + 0.07*sqrt(We),-6);


  /*---- Extended drag coefficient ----*/
  if ( We <= 12 )
    Cd = (1.0 - eccentricity)*Cd_sphere + eccentricity*Cd_disk;
  else if ( We > 12 )
    Cd = Cd_disk;

//  Cd = Cd_sphere;
  Cd = (Re) ? 24/Re * (1 + 0.15*pow(Re,0.687) + 0.0175/(1+ 4.25e4/pow(Re,1.16))) : 1.0;
//  Cd = 24/Re;


  /*---- Compute residual  ----*/

  su2double tau = (Re) ? (4*rho*d*d) / (3*mu*Re*Cd) : 1.0;

  return tau;

}

void CPTSolver::PrintVerificationError(const CConfig *config) const {

  if ((rank != MASTER_NODE) || (MGLevel != MESH_0)) return;

  if (config && !config->GetDiscrete_Adjoint()) {

    cout.precision(5);
    cout.setf(ios::scientific, ios::floatfield);

    cout << endl   << "------------------------ Global Error Analysis --------------------------" << endl;

    cout << setw(20) << "RMS Error  [alpha]: " << setw(12) << VerificationSolution->GetError_RMS(0) << "     | ";
    cout << setw(20) << "Max Error  [alpha]: " << setw(12) << VerificationSolution->GetError_Max(0);
    cout << endl;

    cout << setw(20) << "RMS Error [alphaU]: " << setw(12) << VerificationSolution->GetError_RMS(1) << "     | ";
    cout << setw(20) << "Max Error [alphaU]: " << setw(12) << VerificationSolution->GetError_Max(1);
    cout << endl;

    cout << setw(20) << "RMS Error [alphaV]: " << setw(12) << VerificationSolution->GetError_RMS(2) << "     | ";
    cout << setw(20) << "Max Error [alphaV]: " << setw(12) << VerificationSolution->GetError_Max(2);
    cout << endl;

    if (nDim == 3) {
      cout << setw(20) << "RMS Error [alphaW]: " << setw(12) << VerificationSolution->GetError_RMS(3) << "     | ";
      cout << setw(20) << "Max Error [alphaW]: " << setw(12) << VerificationSolution->GetError_Max(3);
      cout << endl;
    }

    cout << "-------------------------------------------------------------------------" << endl << endl;
    cout.unsetf(ios_base::floatfield);
  }
}


void CPTSolver::ComputeVerificationError(CGeometry* geometry, CConfig* config) {

  /*--- The errors only need to be computed on the finest grid. ---*/
  if (MGLevel != MESH_0) return;

  /*--- If this is a verification case, we can compute the global
   error metrics by using the difference between the local error
   and the known solution at each DOF. This is then collected into
   RMS (L2) and maximum (Linf) global error norms. From these
   global measures, one can compute the order of accuracy. ---*/

  bool write_heads =
      ((((config->GetInnerIter() % (config->GetWrt_Con_Freq() * 40)) == 0) && (config->GetInnerIter() != 0)) ||
       (config->GetInnerIter() == 1));
  if (!write_heads) return;

  /*--- Check if there actually is an exact solution for this
        verification case, if computed at all. ---*/
  if (VerificationSolution && VerificationSolution->ExactSolutionKnown()) {
    /*--- Get the physical time if necessary. ---*/
    su2double time = 0.0;
    if (config->GetTime_Marching()) time = config->GetPhysicalTime();

    /*--- Reset the global error measures to zero. ---*/
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      VerificationSolution->SetError_RMS(iVar, 0.0);
      VerificationSolution->SetError_Max(iVar, 0.0, 0);
    }

    /*--- Loop over all owned points. ---*/
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /* Set the pointers to the coordinates and solution of this DOF. */
      const su2double* coor = geometry->nodes->GetCoord(iPoint);
      su2double* solDOF = nodes->GetSolution(iPoint);

      /* Get local error from the verification solution class. */
      vector<su2double> error(nVar, 0.0);
      VerificationSolution->GetLocalError(coor, time, solDOF, error.data());

      /* Increment the global error measures */
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        VerificationSolution->AddError_RMS(iVar, error[iVar] * error[iVar]);
//        VerificationSolution->AddError_RMS(iVar, error[iVar] * error[iVar] * geometry->nodes->GetVolume(iPoint));
        VerificationSolution->AddError_Max(iVar, fabs(error[iVar]), geometry->nodes->GetGlobalIndex(iPoint),
                                           geometry->nodes->GetCoord(iPoint));
      }
    }

    /* Finalize the calculation of the global error measures. */
    VerificationSolution->SetVerificationError(geometry->GetGlobal_nPointDomain(), config);

    /*--- Screen output of the error metrics. This can be improved
     once the new output classes are in place. ---*/

    PrintVerificationError(config);
  }

}
