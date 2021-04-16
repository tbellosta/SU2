/*!
 * \file CPTOutput.cpp
 * \brief Main subroutines for the particle tracking output
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


#include "../../include/output/CPTOutput.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"
#include "../../include/solvers/CPTSolver.hpp"

CPTOutput::CPTOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

  multiZone = config->GetMultizone_Problem();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_ALPHA");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_ALPHA");
    nRequestedScreenFields = requestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Particle Tracking)";
  multiZoneHeaderString = ss.str();



  
  /*--- Set the volume filename --- */

  volumeFilename = config->GetVolume_FileName_PT();
  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfCoeff_FileName_PT();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_FileName_PT();

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_ALPHA");

  minLogResidual = config->GetMinLogResidual_PT();
  


}

CPTOutput::~CPTOutput(void) = default;

void CPTOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  CSolver* PT_solver  = solver[PT_SOL];
  

  SetHistoryOutputValue("RMS_ALPHA", log10(PT_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_ALPHA-U", log10(PT_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_ALPHA-V", log10(PT_solver->GetRes_RMS(2)));
  if (nDim == 3) SetHistoryOutputValue("RMS_ALPHA-W", log10(PT_solver->GetRes_RMS(3)));

  SetHistoryOutputValue("MAX_ALPHA", log10(PT_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_ALPHA-U", log10(PT_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_ALPHA-V", log10(PT_solver->GetRes_Max(2)));
  if (nDim == 3) SetHistoryOutputValue("MAX_ALPHA-W", log10(PT_solver->GetRes_Max(3)));


  SetHistoryOutputValue("LINSOL_ITER", PT_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(PT_solver->GetResLinSolver()));
  SetHistoryOutputValue("CFL_NUMBER", config->GetCFL(MESH_0));

}


void CPTOutput::SetHistoryOutputFields(CConfig *config){

  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");

  AddHistoryOutput("RMS_ALPHA", "rms[alpha]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the volume fraction", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("RMS_ALPHA-U", "rms[alpha*u]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the momentum", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("RMS_ALPHA-V", "rms[alpha*v]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the momentum", HistoryFieldType::RESIDUAL);
  if (nDim == 3) AddHistoryOutput("RMS_ALPHA-W", "rms[alpha*w]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the momentum", HistoryFieldType::RESIDUAL);

  AddHistoryOutput("MAX_ALPHA", "max[alpha]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the volume fraction", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("MAX_ALPHA-U", "max[alpha*u]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the momentum", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("MAX_ALPHA-V", "max[alpha*v]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the momentum", HistoryFieldType::RESIDUAL);
  if (nDim == 3) AddHistoryOutput("MAX_ALPHA-W", "max[alpha*w]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the momentum", HistoryFieldType::RESIDUAL);

  AddHistoryOutput("CFL_NUMBER", "CFL number", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current value of the CFL number");

}


void CPTOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES","z-component of the coordinate vector");

  // SOLUTION
  AddVolumeOutput("ALPHA", "alpha", "SOLUTION", "Volume fraction");
  AddVolumeOutput("ALPHA-U", "alphaU", "SOLUTION", "Momentum-X");
  AddVolumeOutput("ALPHA-V", "alphaV", "SOLUTION", "Momentum-Y");
  if (nDim == 3) AddVolumeOutput("ALPHA-W", "alphaW", "SOLUTION", "Momentum-Z");

  // PRIMITIVES
  AddVolumeOutput("PART-U", "Up", "PRIMITIVE", "Velocity-X");
  AddVolumeOutput("PART-V", "Vp", "PRIMITIVE", "Velocity-Y");
  AddVolumeOutput("PART-RE", "Re", "PRIMITIVE", "Reynolds-P");
  if (nDim == 3) AddVolumeOutput("PART-W", "Wp", "PRIMITIVE", "Velocity-Z");

  // Residuals
  AddVolumeOutput("RES_ALPHA", "Residual_alpha", "RESIDUAL", "Residual of the volume fraction");
  AddVolumeOutput("RES_ALPHA-U", "Residual_alphaU", "RESIDUAL", "Residual of the MOMENTUM-X");
  AddVolumeOutput("RES_ALPHA-V", "Residual_alphaV", "RESIDUAL", "Residual of the MOMENTUM-Y");
  if (nDim == 3) AddVolumeOutput("RES_ALPHA-W", "Residual_alphaW", "RESIDUAL", "Residual of the MOMENTUM-Z");

  // Mesh quality metrics, computed in CPhysicalGeometry::ComputeMeshQualityStatistics.
  AddVolumeOutput("ORTHOGONALITY", "Orthogonality", "MESH_QUALITY", "Orthogonality Angle (deg.)");
  AddVolumeOutput("ASPECT_RATIO",  "Aspect_Ratio",  "MESH_QUALITY", "CV Face Area Aspect Ratio");
  AddVolumeOutput("VOLUME_RATIO",  "Volume_Ratio",  "MESH_QUALITY", "CV Sub-Volume Ratio");

  AddVolumeOutput("BETA", "CollEff", "SOLUTION", "Collection Efficiency");
  AddVolumeOutput("BETA_CORRECTED", "CollEff_Corrected", "SOLUTION", "Collection Efficiency Corrected for Splashing");
  
  if(config->GetMultiBin()){
    AddVolumeOutput("BETA_TOT", "CollEff_tot", "SOLUTION", "Total Collection Efficiency");
    AddVolumeOutput("BETA_CORRECTED_TOT", "CollEff_Corrected_tot", "SOLUTION", " Total Collection Efficiency Corrected for Splashing");
  }

  if (config->GetKind_ConvNumScheme_PT() == SPACE_CENTERED)
  AddVolumeOutput("SENSOR", "Sensor", "SOLUTION", "Sensor");


  // LIMITERS
  if (config->GetKind_SlopeLimit_PT() != NO_LIMITER && config->GetMUSCL_PT()) {
    AddVolumeOutput("LIM-U", "limU", "LIMITER", "Limiter-U");
    AddVolumeOutput("LIM-V", "limV", "LIMITER", "Limiter-V");
    AddVolumeOutput("LIM-A", "limA", "LIMITER", "Limiter-alpha");
    if (nDim == 3) AddVolumeOutput("LIM-W", "limW", "LIMITERS", "Limiter-W");
  }
}


void CPTOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
  CVariable* Node_PT = nullptr;
  
  Node_PT = solver[PT_SOL]->GetNodes();
  
  CVariable* Node_Flow = solver[FLOW_SOL]->GetNodes();
  CPoint*    Node_Geo  = geometry->nodes;

  // Grid coordinates
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));

  // SOLUTION
  SetVolumeOutputValue("ALPHA", iPoint, Node_PT->GetSolution(iPoint, 0));
  SetVolumeOutputValue("ALPHA-U", iPoint, Node_PT->GetSolution(iPoint, 1));
  SetVolumeOutputValue("ALPHA-V", iPoint, Node_PT->GetSolution(iPoint, 2));
  if (nDim == 3) SetVolumeOutputValue("ALPHA-W", iPoint, Node_PT->GetSolution(iPoint, 3));

  // PRIMITIVES
  SetVolumeOutputValue("PART-U", iPoint, Node_PT->GetPrimitive(iPoint, 1));
  SetVolumeOutputValue("PART-V", iPoint, Node_PT->GetPrimitive(iPoint, 2));
  if (nDim == 3) SetVolumeOutputValue("PART-W", iPoint, Node_PT->GetPrimitive(iPoint, 3));

  SetVolumeOutputValue("RES_ALPHA", iPoint, solver[PT_SOL]->LinSysRes(iPoint, 0));
  SetVolumeOutputValue("RES_ALPHA-U", iPoint, solver[PT_SOL]->LinSysRes(iPoint, 1));
  SetVolumeOutputValue("RES_ALPHA-V", iPoint, solver[PT_SOL]->LinSysRes(iPoint, 2));
  //  SetVolumeOutputValue("RES_ALPHA-U", iPoint, Node_PT->GetGradient_Primitive(iPoint,1,0));
  //  SetVolumeOutputValue("RES_ALPHA-V", iPoint, Node_PT->GetGradient_Primitive(iPoint,1,1));
  if (nDim == 3) SetVolumeOutputValue("RES_ALPHA-W", iPoint, solver[PT_SOL]->LinSysRes(iPoint, 3));
  

  // Mesh quality metrics
  if (config->GetWrt_MeshQuality()) {
    SetVolumeOutputValue("ORTHOGONALITY", iPoint, geometry->Orthogonality[iPoint]);
    SetVolumeOutputValue("ASPECT_RATIO",  iPoint, geometry->Aspect_Ratio[iPoint]);
    SetVolumeOutputValue("VOLUME_RATIO",  iPoint, geometry->Volume_Ratio[iPoint]);
  }

  const su2double rho = 1000;
  const su2double sigma = 0.0756;
//  const su2double d = config->GetParticle_Size();
  const su2double d = 2e-5;
  su2double mu = 18.03e-6;

  su2double *uFlow = &Node_Flow->GetPrimitive(iPoint)[1], rhoFlow = Node_Flow->GetPrimitive(iPoint)[nDim+2];

  su2double V_mag = 0.0;
  for (int iDim = 0; iDim < nDim; ++iDim) {
    V_mag += pow(uFlow[iDim]-Node_PT->GetPrimitive(iPoint,iDim+1),2);
  }

  su2double Re = rhoFlow*d*V_mag / mu;

  SetVolumeOutputValue("PART-RE", iPoint, Re);

  if (config->GetKind_ConvNumScheme_PT() == SPACE_CENTERED)
    SetVolumeOutputValue("SENSOR",  iPoint, Node_PT->GetSensor(iPoint));

  if (config->GetKind_SlopeLimit_PT() != NO_LIMITER && config->GetMUSCL_PT()) {
    SetVolumeOutputValue("LIM-U", iPoint, Node_PT->GetLimiter_Primitive(iPoint,1));
    SetVolumeOutputValue("LIM-V", iPoint, Node_PT->GetLimiter_Primitive(iPoint,2));
    SetVolumeOutputValue("LIM-A", iPoint, Node_PT->GetLimiter_Primitive(iPoint,0));
    if (nDim == 3) SetVolumeOutputValue("LIM-W", iPoint, Node_PT->GetLimiter_Primitive(iPoint,3));
  }


}

void CPTOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){
  
  SetVolumeOutputValue("BETA", iPoint, fmax(0,solver[PT_SOL]->GetCollectionEfficiency(iMarker, iVertex)));
    
  CPTSolver *PTSolver = dynamic_cast<CPTSolver*>(solver[PT_SOL]);
  if(config->GetSplashingPT() && PTSolver->CollectionEfficiencyCorrectedSplashing !=nullptr){

    SetVolumeOutputValue("BETA_CORRECTED", iPoint, fmax(0,PTSolver->CollectionEfficiencyCorrectedSplashing[iMarker][iVertex]));
  }
  if(config->GetMultiBin() && PTSolver->CollectionEfficiencyTOT !=nullptr){
    SetVolumeOutputValue("BETA_TOT", iPoint, fmax(0,PTSolver->CollectionEfficiencyTOT[iMarker][iVertex]));
    if(config->GetSplashingPT() && PTSolver->CollectionEfficiencyCorrectedSplashingTOT !=nullptr){

      SetVolumeOutputValue("BETA_CORRECTED_TOT", iPoint, fmax(0,PTSolver->CollectionEfficiencyCorrectedSplashingTOT[iMarker][iVertex]));
    }
  }
  
}

