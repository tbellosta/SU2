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

CPTOutput::CPTOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

  multiZone = config->GetMultizone_Problem();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_RES");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_DENSITY");
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

  volumeFilename = config->GetVolume_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_FileName();

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_DENSITY");


}

CPTOutput::~CPTOutput(void) = default;

void CPTOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  CSolver* PT_solver = solver[PT_SOL];

  SetHistoryOutputValue("RMS_DENSITY", log10(PT_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_MOMENTUM-X", log10(PT_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_MOMENTUM-Y", log10(PT_solver->GetRes_RMS(2)));
  if (nDim == 3) SetHistoryOutputValue("RMS_MOMENTUM-Z", log10(PT_solver->GetRes_RMS(3)));

  SetHistoryOutputValue("MAX_DENSITY", log10(PT_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_MOMENTUM-X", log10(PT_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_MOMENTUM-Y", log10(PT_solver->GetRes_Max(2)));
  if (nDim == 3) SetHistoryOutputValue("MAX_MOMENTUM-Z", log10(PT_solver->GetRes_Max(3)));


  SetHistoryOutputValue("LINSOL_ITER", PT_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(PT_solver->GetResLinSolver()));
  SetHistoryOutputValue("CFL_NUMBER", config->GetCFL(MESH_0));

}


void CPTOutput::SetHistoryOutputFields(CConfig *config){

  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");

  AddHistoryOutput("RMS_DENSITY", "rms[alpha]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the volume fraction", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("RMS_MOMENTUM-X", "rms[alpha*u]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the momentum", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("RMS_MOMENTUM-Y", "rms[alpha*v]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the momentum", HistoryFieldType::RESIDUAL);
  if (nDim == 3) AddHistoryOutput("RMS_MOMENTUM-Z", "rms[alpha*w]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the momentum", HistoryFieldType::RESIDUAL);

  AddHistoryOutput("MAX_DENSITY", "max[alpha]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the volume fraction", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("MAX_MOMENTUM-X", "max[alpha*u]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the momentum", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("MAX_MOMENTUM-Y", "max[alpha*v]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the momentum", HistoryFieldType::RESIDUAL);
  if (nDim == 3) AddHistoryOutput("MAX_MOMENTUM-Z", "max[alpha*w]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the momentum", HistoryFieldType::RESIDUAL);

  AddHistoryOutput("CFL_NUMBER", "CFL number", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current value of the CFL number");

}


void CPTOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES","z-component of the coordinate vector");

  // SOLUTION
  AddVolumeOutput("DENSITY", "alpha", "SOLUTION", "Volume fraction");
  AddVolumeOutput("MOMENTUM-X", "alphaU", "SOLUTION", "Momentum-X");
  AddVolumeOutput("MOMENTUM-Y", "alphaV", "SOLUTION", "Momentum-Y");
  if (nDim == 3) AddVolumeOutput("MOMENTUM-Z", "alphaW", "SOLUTION", "Momentum-Z");

//  // Primitives
//  AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heatflux");

  // Residuals
  AddVolumeOutput("RES_DENSITY", "Residual_alpha", "RESIDUAL", "Residual of the volume fraction");
  AddVolumeOutput("RES_MOMENTUM-X", "Residual_alphaU", "RESIDUAL", "Residual of the MOMENTUM-X");
  AddVolumeOutput("RES_MOMENTUM-Y", "Residual_alphaV", "RESIDUAL", "Residual of the MOMENTUM-Y");
  if (nDim == 3) AddVolumeOutput("RES_MOMENTUM-Z", "Residual_alphaW", "RESIDUAL", "Residual of the MOMENTUM-Z");

  // Mesh quality metrics, computed in CPhysicalGeometry::ComputeMeshQualityStatistics.
  AddVolumeOutput("ORTHOGONALITY", "Orthogonality", "MESH_QUALITY", "Orthogonality Angle (deg.)");
  AddVolumeOutput("ASPECT_RATIO",  "Aspect_Ratio",  "MESH_QUALITY", "CV Face Area Aspect Ratio");
  AddVolumeOutput("VOLUME_RATIO",  "Volume_Ratio",  "MESH_QUALITY", "CV Sub-Volume Ratio");

}


void CPTOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_PT = solver[PT_SOL]->GetNodes();
  CPoint*    Node_Geo  = geometry->nodes;

  // Grid coordinates
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));

  // SOLUTION
  SetVolumeOutputValue("DENSITY", iPoint, Node_PT->GetSolution(iPoint, 0));
  SetVolumeOutputValue("MOMENTUM-X", iPoint, Node_PT->GetSolution(iPoint, 1));
  SetVolumeOutputValue("MOMENTUM-Y", iPoint, Node_PT->GetSolution(iPoint, 2));
  if (nDim == 3) SetVolumeOutputValue("MOMENTUM-Z", iPoint, Node_PT->GetSolution(iPoint, 3));

  // Residuals
  SetVolumeOutputValue("RES_DENSITY", iPoint, solver[PT_SOL]->LinSysRes(iPoint, 0));
  SetVolumeOutputValue("RES_MOMENTUM-X", iPoint, solver[PT_SOL]->LinSysRes(iPoint, 1));
  SetVolumeOutputValue("RES_MOMENTUM-Y", iPoint, solver[PT_SOL]->LinSysRes(iPoint, 2));
  if (nDim == 3) SetVolumeOutputValue("RES_MOMENTUM-Z", iPoint, solver[PT_SOL]->LinSysRes(iPoint, 3));

  // Mesh quality metrics
  if (config->GetWrt_MeshQuality()) {
    SetVolumeOutputValue("ORTHOGONALITY", iPoint, geometry->Orthogonality[iPoint]);
    SetVolumeOutputValue("ASPECT_RATIO",  iPoint, geometry->Aspect_Ratio[iPoint]);
    SetVolumeOutputValue("VOLUME_RATIO",  iPoint, geometry->Volume_Ratio[iPoint]);
  }

}

void CPTOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){

  /** @todo Here add collection efficiency output **/
//  SetVolumeOutputValue("BETA", iPoint, solver[PT_SOL]->GetCollectionEfficiency(iMarker, iVertex));

}

