/*!
 * \file heat.cpp
 * \brief Implementation of numerics classes for particle tracking.
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

#include "../../include/numerics/particleTracking.hpp"

CConv_PT::CConv_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  /** A grid is defined as dynamic if there's rigid grid movement or
   *  grid deformation AND the problem is time domain **/
  dynamic_grid = config->GetDynamic_Grid();

  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];

  Proj_Flux_i = new su2double [nVar];
  Proj_Flux_j = new su2double [nVar];

  for (int iVar = 0; iVar < nVar; ++iVar) {
    Proj_Jac_i = new su2double*[nVar];
    Proj_Jac_j = new su2double*[nVar];
    for (int jVar = 0; jVar < nVar; ++jVar) {
      Proj_Jac_i[iVar] = new su2double[nVar];
      Proj_Jac_j[iVar] = new su2double[nVar];
    }
  }

  Laminar_Viscosity_i = config->GetViscosity_FreeStreamND();
  Laminar_Viscosity_j = config->GetViscosity_FreeStreamND();
}

CConv_PT::~CConv_PT() {

  delete [] Velocity_i;
  delete [] Velocity_j;

  delete [] Proj_Flux_i;
  delete [] Proj_Flux_j;

  for (int iVar = 0; iVar < nVar; ++iVar) {
    delete [] Proj_Jac_i[iVar];
    delete [] Proj_Jac_j[iVar];
  }

  delete [] Proj_Jac_i;
  delete [] Proj_Jac_j;
}


void CConv_PT::GetProjFluxPT(const su2double* VolFraction, const su2double* Vel, const su2double* Norm,
                             su2double* ProjFlux) const {

  su2double au = (*VolFraction) * Vel[0];
  su2double av = (*VolFraction) * Vel[1];

  ProjFlux[0] = au * Norm[0];
  ProjFlux[1] = au * Vel[0] * Norm[0];
  ProjFlux[2] = au * Vel[1] * Norm[0];

  ProjFlux[0] += av * Norm[1];
  ProjFlux[1] += av * Vel[0] * Norm[1];
  ProjFlux[2] += av * Vel[1] * Norm[1];

  if (nDim == 3) {

    su2double aw = (*VolFraction) * Vel[2];

    ProjFlux[3] = au * Vel[2] * Norm[0];
    ProjFlux[3] += av * Vel[2] * Norm[1];

    ProjFlux[0] += aw * Norm[1];
    ProjFlux[1] += aw * Vel[0] * Norm[2];
    ProjFlux[2] += aw * Vel[1] * Norm[2];
    ProjFlux[3] += aw * Vel[2] * Norm[2];

  }

}
void CConv_PT::GetProjFluxJacobianPT(const su2double* VolFraction, const su2double* Vel, const su2double* Norm,
                                     su2double** ProjJac) const {

  su2double u = Vel[0];
  su2double v = Vel[1];

  ProjJac[0][0]  =    0 * Norm[0];   ProjJac[0][1]  =   1 * Norm[0];  ProjJac[0][2]  =   0 * Norm[0];
  ProjJac[1][0]  = -u*u * Norm[0];   ProjJac[1][1]  = 2*u * Norm[0];  ProjJac[1][2]  =   0 * Norm[0];
  ProjJac[2][0]  = -u*v * Norm[0];   ProjJac[2][1]  =   v * Norm[0];  ProjJac[2][2]  =   u * Norm[0];

  ProjJac[0][0] +=    0 * Norm[1];   ProjJac[0][1] +=   0 * Norm[1];  ProjJac[0][2] +=   1 * Norm[1];
  ProjJac[1][0] += -u*v * Norm[1];   ProjJac[1][1] +=   v * Norm[1];  ProjJac[1][2] +=   u * Norm[1];
  ProjJac[2][0] += -v*v * Norm[1];   ProjJac[2][1] +=   0 * Norm[1];  ProjJac[2][2] += 2*v * Norm[1];


  if (nDim == 3) {
    su2double w = Vel[2];

                                                                                                       ProjJac[0][3]  =   0 * Norm[0];
                                                                                                       ProjJac[1][3]  =   0 * Norm[0];
                                                                                                       ProjJac[2][3]  =   0 * Norm[0];
    ProjJac[3][0]  = -u*w * Norm[0];   ProjJac[3][1]  = w * Norm[0];   ProjJac[3][2] = 0 * Norm[0];    ProjJac[3][3]  =   u * Norm[0];

                                                                                                       ProjJac[0][3] +=   0 * Norm[1];
                                                                                                       ProjJac[1][3] +=   0 * Norm[1];
                                                                                                       ProjJac[2][3] +=   0 * Norm[1];
    ProjJac[3][0] += -v*w * Norm[1];   ProjJac[3][1] += 0 * Norm[1];   ProjJac[3][2] += w * Norm[1];   ProjJac[3][3] +=   v * Norm[1];

    ProjJac[0][0] +=    0 * Norm[2];   ProjJac[0][1] += 0 * Norm[2];   ProjJac[0][2] += 0 * Norm[2];   ProjJac[0][3] +=   1 * Norm[2];
    ProjJac[1][0] += -u*w * Norm[2];   ProjJac[1][1] += w * Norm[2];   ProjJac[1][2] += 0 * Norm[2];   ProjJac[1][3] +=   u * Norm[2];
    ProjJac[2][0] += -v*w * Norm[2];   ProjJac[2][1] += 0 * Norm[2];   ProjJac[2][2] += w * Norm[2];   ProjJac[2][3] +=   v * Norm[2];
    ProjJac[3][0] += -w*w * Norm[2];   ProjJac[3][1] += 0 * Norm[2];   ProjJac[3][2] += 0 * Norm[2];   ProjJac[3][3] += 2*w * Norm[2];
  }


}






CUpwRusanov_PT::CUpwRusanov_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
    CConv_PT(val_nDim, val_nVar, config) {}


void CUpwRusanov_PT::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                   su2double **val_Jacobian_j, CConfig *config) {


  su2double projVel_i = 0;
  su2double projVel_j = 0;
  su2double lambda;

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    projVel_i += Velocity_i[iDim]*Normal[iDim];
    projVel_j += Velocity_j[iDim]*Normal[iDim];
  }

  lambda = fmax(fabs(projVel_i),fabs(projVel_j));

  GetProjFluxPT(&Density_i, Velocity_i, Normal, Proj_Flux_i);
  GetProjFluxPT(&Density_j, Velocity_j, Normal, Proj_Flux_j);

  for (int iVar = 0; iVar < nVar; ++iVar) {
    val_residual[iVar] = 0.5*(Proj_Flux_i[iVar]+Proj_Flux_j[iVar]) - 0.5*lambda*(U_j[iVar] - U_i[iVar]);
  }


  if (implicit) {

    GetProjFluxJacobianPT(&Density_i, Velocity_i, Normal, Proj_Jac_i);
    GetProjFluxJacobianPT(&Density_j, Velocity_j, Normal, Proj_Jac_j);

    for (int iVar = 0; iVar < nVar; ++iVar) {
      for (int jVar = 0; jVar < nVar; ++jVar) {
        val_Jacobian_i[iVar][jVar] = 0.5*Proj_Jac_i[iVar][jVar];
        val_Jacobian_j[iVar][jVar] = 0.5*Proj_Jac_j[iVar][jVar];
      }
    }

    for (int iVar = 0; iVar < nVar; ++iVar) {
      val_Jacobian_i[iVar][iVar] += 0.5*lambda;
      val_Jacobian_j[iVar][iVar] -= 0.5*lambda;
    }

  }


}