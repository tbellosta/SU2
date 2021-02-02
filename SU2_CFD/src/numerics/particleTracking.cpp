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

  implicit = (config->GetKind_TimeIntScheme_PT() == EULER_IMPLICIT);
  /** A grid is defined as dynamic if there's rigid grid movement or
   *  grid deformation AND the problem is time domain **/
  dynamic_grid = config->GetDynamic_Grid();

  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];

  Conservatives_i = new su2double [nVar];
  Conservatives_j = new su2double [nVar];


  Proj_Flux_i = new su2double [nVar];
  Proj_Flux_j = new su2double [nVar];

  Proj_Jac_i = new su2double*[nVar];
  Proj_Jac_j = new su2double*[nVar];
  for (int jVar = 0; jVar < nVar; ++jVar) {
    Proj_Jac_i[jVar] = new su2double[nVar];
    Proj_Jac_j[jVar] = new su2double[nVar];
  }

  Laminar_Viscosity_i = config->GetViscosity_FreeStreamND();
  Laminar_Viscosity_j = config->GetViscosity_FreeStreamND();
}

CConv_PT::~CConv_PT() {

  delete [] Velocity_i;
  delete [] Velocity_j;

  delete [] Conservatives_j;
  delete [] Conservatives_i;

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

    ProjFlux[0] += aw * Norm[2];
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
    CConv_PT(val_nDim, val_nVar, config) { }


void CUpwRusanov_PT::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                   su2double **val_Jacobian_j, CConfig *config) {
#define SMALL 1e-50

  su2double projVel_i = 0;
  su2double projVel_j = 0;
  su2double lambda;

  Density_i = V_i[0];
  Density_j = V_j[0];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    projVel_i += Velocity_i[iDim]*Normal[iDim];
    projVel_j += Velocity_j[iDim]*Normal[iDim];
  }

  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;

  for (iDim = 0; iDim < nDim; ++iDim) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }

  lambda = fmax(fabs(projVel_i),fabs(projVel_j));

  GetProjFluxPT(&Density_i, Velocity_i, Normal, Proj_Flux_i);
  GetProjFluxPT(&Density_j, Velocity_j, Normal, Proj_Flux_j);

  for (int iVar = 0; iVar < nVar; ++iVar) {
    val_residual[iVar] = 0.5*(Proj_Flux_i[iVar]+Proj_Flux_j[iVar]) - 0.5*lambda*(Conservatives_j[iVar] - Conservatives_i[iVar]);
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


CUpwGodunov_PT::CUpwGodunov_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
    CConv_PT(val_nDim, val_nVar, config) { }


void CUpwGodunov_PT::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                     su2double **val_Jacobian_j, CConfig *config) {
#define SMALL 1e-50

  su2double projVel_i = 0;
  su2double projVel_j = 0;
  su2double uDelta;

  Density_i = V_i[0];
  Density_j = V_j[0];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    projVel_i += Velocity_i[iDim]*Normal[iDim];
    projVel_j += Velocity_j[iDim]*Normal[iDim];
  }

  uDelta = (sqrt(Density_i)*projVel_i + sqrt(Density_j)*projVel_j) / (sqrt(Density_i) + sqrt(Density_j) + EPS);

  GetProjFluxPT(&Density_i, Velocity_i, Normal, Proj_Flux_i);
  GetProjFluxPT(&Density_j, Velocity_j, Normal, Proj_Flux_j);

  if (projVel_i >= projVel_j) {
    if (uDelta > 0) {
      for (int iVar = 0; iVar < nVar; ++iVar)
        val_residual[iVar] = Proj_Flux_i[iVar];
    } else if (uDelta < 0) {
      for (int iVar = 0; iVar < nVar; ++iVar)
        val_residual[iVar] = Proj_Flux_j[iVar];
    } else {
      for (int iDim = 0; iDim < nDim; ++iDim) {
        val_residual[iDim+1] = Proj_Flux_i[iDim+1];
      }
      val_residual[0] = 0.5 * (Proj_Flux_j[0] + Proj_Flux_i[0]);
    }
  } else {
    if (projVel_i > 0){
      for (int iVar = 0; iVar < nVar; ++iVar)
        val_residual[iVar] = Proj_Flux_i[iVar];
    } else if (projVel_j < 0){
      for (int iVar = 0; iVar < nVar; ++iVar)
        val_residual[iVar] = Proj_Flux_j[iVar];
    } else {
      for (int iVar = 0; iVar < nVar; ++iVar)
        val_residual[iVar] = 0.0;
    }
  }


  if (implicit) {

    GetProjFluxJacobianPT(&Density_i, Velocity_i, Normal, Proj_Jac_i);
    GetProjFluxJacobianPT(&Density_j, Velocity_j, Normal, Proj_Jac_j);

    for (int iVar = 0; iVar < nVar; ++iVar) {
      for (int jVar = 0; jVar < nVar; ++jVar) {
        val_Jacobian_i[iVar][jVar] = 0.0;
        val_Jacobian_j[iVar][jVar] = 0.0;
      }
    }



    if (projVel_i >= projVel_j) {
      if (uDelta > 0) {
        for (int iVar = 0; iVar < nVar; ++iVar)
          for (int jVar = 0; jVar < nVar; ++jVar)
            val_Jacobian_i[iVar][jVar] = Proj_Jac_i[iVar][jVar];

      } else if (uDelta < 0) {
        for (int iVar = 0; iVar < nVar; ++iVar)
          for (int jVar = 0; jVar < nVar; ++jVar)
            val_Jacobian_j[iVar][jVar] = Proj_Jac_j[iVar][jVar];

      } else {
        for (int iVar = 1; iVar < nVar; ++iVar)
          for (int jVar = 0; jVar < nVar; ++jVar)
            val_Jacobian_i[iVar][jVar] = Proj_Jac_i[iVar][jVar];


        for (int jVar = 0; jVar < nVar; ++jVar) {
          val_Jacobian_i[0][jVar] = 0.5*Proj_Jac_i[0][jVar];
          val_Jacobian_j[0][jVar] = 0.5*Proj_Jac_j[0][jVar];
        }

      }
    } else {
      if (projVel_i > 0){
        for (int iVar = 0; iVar < nVar; ++iVar)
          for (int jVar = 0; jVar < nVar; ++jVar)
            val_Jacobian_i[iVar][jVar] = Proj_Jac_i[iVar][jVar];
      } else if (projVel_j < 0){
        for (int iVar = 0; iVar < nVar; ++iVar)
          for (int jVar = 0; jVar < nVar; ++jVar)
            val_Jacobian_j[iVar][jVar] = Proj_Jac_j[iVar][jVar];
      }
    }

  }


}


CUpwFDS_PT::CUpwFDS_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
    CConv_PT(val_nDim, val_nVar, config) { }

void CUpwFDS_PT::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  su2double projVel_i = 0;
  su2double projVel_j = 0;
  su2double uRoe;

  Density_i = V_i[0];
  Density_j = V_j[0];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    projVel_i += Velocity_i[iDim]*Normal[iDim];
    projVel_j += Velocity_j[iDim]*Normal[iDim];
  }

  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;

  for (iDim = 0; iDim < nDim; ++iDim) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }

  uRoe = (sqrt(Density_i)*projVel_i + sqrt(Density_j)*projVel_j) / (sqrt(Density_i) + sqrt(Density_j));
  uRoe = fabs(uRoe);

//  su2double eps = 1e-2;
//
//  if (uRoe < eps) uRoe = 0.5*(uRoe*uRoe/eps + eps);

  GetProjFluxPT(&Density_i, Velocity_i, Normal, Proj_Flux_i);
  GetProjFluxPT(&Density_j, Velocity_j, Normal, Proj_Flux_j);

  for (int iVar = 0; iVar < nVar; ++iVar) {
    val_residual[iVar] = 0.5 * ((Proj_Flux_i[iVar] + Proj_Flux_j[iVar]) - uRoe * (Conservatives_j[iVar] - Conservatives_i[iVar]));
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

    /*--- Diagonal terms ---*/
    for (int iVar = 0; iVar < nVar; ++iVar) {
      val_Jacobian_i[iVar][iVar] += 0.5 * uRoe;
      val_Jacobian_j[iVar][iVar] -= 0.5 * uRoe;
    }
  }
}



void CSourceDrag::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j,
                                  CConfig* config) {

  su2double *uFlow = &V_i[1];

  su2double Uf[2] = {10, 10};

  uFlow = Uf;

  val_residual[0] = 0.0;

    for (int iDim = 0; iDim < nDim; ++iDim) {
      val_residual[iDim+1] = Volume*U_i[0]*(uFlow[iDim]-U_i[iDim+1]) / ParticleRelaxationTime;
    }


  for (int iVar = 0; iVar < nVar; ++iVar)
    for (int jVar = 0; jVar < nVar; ++jVar)
      val_Jacobian_i[iVar][jVar] = 0.0;



  su2double k = Volume / ParticleRelaxationTime;
  for (int iDim = 0; iDim < nDim; ++iDim) {
    val_Jacobian_i[iDim+1][iDim+1] = -k;
    val_Jacobian_i[iDim+1][0] = uFlow[iDim] * k;
  }


}
CSourceDrag::CSourceDrag(unsigned short val_nDim, unsigned short val_nVar, CConfig* config) :
                        CNumerics(val_nDim, val_nVar, config) {}



CUpwStegWarm_PT::CUpwStegWarm_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig* config)
    : CConv_PT(val_nDim, val_nVar, config) {}


void CUpwStegWarm_PT::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j,
                                      CConfig* config) {
  su2double projVel_i = 0;
  su2double projVel_j = 0;

  Density_i = V_i[0];
  Density_j = V_j[0];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    projVel_i += Velocity_i[iDim]*Normal[iDim];
    projVel_j += Velocity_j[iDim]*Normal[iDim];
  }

  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;

  for (iDim = 0; iDim < nDim; ++iDim) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }

  for (int iVar = 0; iVar < nVar; ++iVar) {
    val_residual[iVar] = 0.5*((projVel_i+fabs(projVel_i))*Conservatives_i[iVar] + (projVel_j-fabs(projVel_j))*Conservatives_j[iVar]);
  }


  if (implicit) {

    for (int iVar = 0; iVar < nVar; ++iVar) {
      for (int jVar = 0; jVar < nVar; ++jVar) {
        val_Jacobian_i[iVar][jVar] = 0;
        val_Jacobian_j[iVar][jVar] = 0;
      }
    }

    for (int iVar = 0; iVar < nVar; ++iVar) {
      val_Jacobian_i[iVar][iVar] = 0.5*(projVel_i+fabs(projVel_i));
      val_Jacobian_j[iVar][iVar] = 0.5*(projVel_j-fabs(projVel_j));
    }

  }

}


CAvgGradCorrected_PT::CAvgGradCorrected_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
    CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_PT() == EULER_IMPLICIT);

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradPrimVar_Edge = new su2double [nVar];
  Proj_Mean_GradPrimVar_Kappa = new su2double [nVar];
  Proj_Mean_GradPrimVar_Corrected = new su2double [nVar];
  Mean_GradPrimVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) Mean_GradPrimVar[iVar] = new su2double [nDim];

}

CAvgGradCorrected_PT::~CAvgGradCorrected_PT(void) {

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradPrimVar_Edge;
  delete [] Proj_Mean_GradPrimVar_Kappa;
  delete [] Proj_Mean_GradPrimVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;

}

void CAvgGradCorrected_PT::ComputeResidual(su2double *val_residual, su2double **Jacobian_i,
                                             su2double **Jacobian_j, CConfig *config) {

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

//  SetNumericalViscosity();
  viscosity_Mean = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);

  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    Proj_Mean_GradPrimVar_Kappa[iVar] = 0.0;

    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPrimVar_Kappa[iVar] += Mean_GradPrimVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradPrimVar_Corrected[iVar] = Proj_Mean_GradPrimVar_Kappa[iVar];
    Proj_Mean_GradPrimVar_Corrected[iVar] -=
        Proj_Mean_GradPrimVar_Edge[iVar]*proj_vector_ij -
                                             (V_j[iVar]-V_i[iVar])*proj_vector_ij;
  }

  val_residual[0] = viscosity_Mean*Proj_Mean_GradPrimVar_Corrected[0];
  for (iDim = 0; iDim < nDim; ++iDim) {
    val_residual[iDim+1] = 0.0;
  }

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    for (iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar) {
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_j[iVar][jVar] = 0.0;
      }

    Jacobian_i[0][0] = -viscosity_Mean * proj_vector_ij;
    Jacobian_j[0][0] = viscosity_Mean * proj_vector_ij;
  }

}

CAvgGrad_PT::CAvgGrad_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig* config)
    : CAvgGradCorrected_PT(val_nDim, val_nVar, config) {}


void CAvgGrad_PT::ComputeResidual(su2double* val_residual, su2double** Jacobian_i, su2double** Jacobian_j,
                                  CConfig* config) {

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  viscosity_Mean = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);

  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPrimVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradPrimVar_Corrected[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_i[iVar][iDim]);
      Proj_Mean_GradPrimVar_Kappa[iVar] += Mean_GradPrimVar[iVar][iDim]*Normal[iDim];
    }
    Proj_Mean_GradPrimVar_Corrected[iVar] = Proj_Mean_GradPrimVar_Kappa[iVar];
  }



  val_residual[0] = viscosity_Mean*Proj_Mean_GradPrimVar_Corrected[0];
  for (iDim = 0; iDim < nDim; ++iDim) {
    val_residual[iDim+1] = 0.0;
  }

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    for (iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar) {
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_j[iVar][jVar] = 0.0;
      }


    Jacobian_i[0][0] = -viscosity_Mean * proj_vector_ij;
    Jacobian_j[0][0] = viscosity_Mean * proj_vector_ij;
  }

}


CCent_PT::CCent_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig* config)
    : CConv_PT(val_nDim, val_nVar, config) {

  MeanVelocity = new su2double [nVar];
  Diff_U       = new su2double [nVar];
  Diff_Lapl    = new su2double [nVar];

  Param_p = 0.3;
  Param_Kappa_2 = 0.5;
  Param_Kappa_4 = 0.02;

  Param_Kappa_2 = 0.1;
  Param_Kappa_4 = 0.01;

  fix_factor = 1.0;
}

CCent_PT::~CCent_PT() {
  delete [] MeanVelocity;
  delete [] Diff_U;
  delete [] Diff_Lapl;
}


void CCent_PT::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j,
                               CConfig* config) {


  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);


  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/

  Density_i  = V_i[0];                       Density_j  = V_j[0];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
  }

  /*--- Recompute conservative variables ---*/

  Conservatives_i[0] = Density_i; Conservatives_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }

  /*--- Compute mean values of the variables ---*/

  MeanDensity = 0.5*(Density_i+Density_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);


  /*--- Residual of the inviscid flux ---*/

  GetProjFluxPT(&Density_i, Velocity_i, Normal, Proj_Flux_i);
  GetProjFluxPT(&Density_j, Velocity_j, Normal, Proj_Flux_j);

  /*--- Jacobians of the inviscid flux, scale = 0.5 because ProjFlux ~ 0.5*(fc_i+fc_j)*Normal ---*/

  if (implicit) {
    GetProjFluxJacobianPT(&Density_i, Velocity_i, Normal, Proj_Jac_i);
    GetProjFluxJacobianPT(&Density_j, Velocity_j, Normal, Proj_Jac_j);
  }

  for (iVar = 0; iVar < nVar; ++iVar) {
    val_residual[iVar] = 0.5*(Proj_Flux_i[iVar] + Proj_Flux_j[iVar]);
    for (jVar = 0; jVar < nVar; ++jVar) {
      val_Jacobian_i[iVar][jVar] = 0.5 * Proj_Jac_i[iVar][jVar];
      val_Jacobian_j[iVar][jVar] = 0.5 * Proj_Jac_j[iVar][jVar];
    }
  }

  /*--- Adjustment due to grid motion ---*/

  su2double ProjGridVel = 0.0;
  if (dynamic_grid) {
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjGridVel * 0.5*(Conservatives_i[iVar] + Conservatives_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjGridVel;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjGridVel;
      }
    }
  }

  return;

  /*--- Compute the local spectral radius and the stretching factor ---*/

  su2double ProjVelocity_i = 0.0, ProjVelocity_j = 0.0;
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  /*--- Adjustment due to mesh motion ---*/

  if (dynamic_grid) {
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }

  /*--- Dissipation term ---*/

  su2double Local_Lambda_i = (fabs(ProjVelocity_i));
  su2double Local_Lambda_j = (fabs(ProjVelocity_j));
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  su2double Phi_i = pow(Lambda_i/(4.0*MeanLambda + EPS), Param_p);
  su2double Phi_j = pow(Lambda_j/(4.0*MeanLambda + EPS), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j + EPS);

  /*--- Compute differences btw. conservative variables, with a correction for enthalpy ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_U[iVar] = Conservatives_i[iVar]-Conservatives_j[iVar];
  }

  DissipationTerm(val_residual, val_Jacobian_i, val_Jacobian_j);

}

void CCent_PT::DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) {


  su2double DT_i, DT_j, DT;
  su2double ProjVel_i = 0.0, ProjVel_j = 0.0, Area = 0.0;

  for (iDim = 0; iDim < nDim; ++iDim) {
    ProjVel_i += Velocity_i[iDim] * Normal[iDim];
    ProjVel_j += Velocity_j[iDim] * Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  /*-- CV Volume is set in the place of laminar viscosity --*/
  DT_i = Laminar_Viscosity_i / ProjVel_i;
  DT_j = Laminar_Viscosity_j / ProjVel_j;
  DT = 0.5*(DT_i + DT_j);

  /*--- Compute differences btw. Laplacians ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
  }

  /*--- Compute dissipation coefficients ---*/

  sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  sc4 = sc2*sc2/4.0;

  sc2 = sc4 = 1.0;

  Epsilon_2 = Param_Kappa_2 * 0.5 * (Sensor_i+Sensor_j) / DT;
  Epsilon_4 = max(0.0, Param_Kappa_4/DT - Epsilon_2);
//  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
//  Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;

  /*--- Compute viscous part of the residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += (Epsilon_2*Diff_U[iVar] - Epsilon_4*Diff_Lapl[iVar])*Area;
//    val_residual[iVar] += (Epsilon_2*Diff_U[iVar] - Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda*10;

  /*--- Jacobian computation ---*/

  if (implicit) {

    cte_0 = (Epsilon_2 + Epsilon_4*su2double(Neighbor_i+1))*Area;
//    cte_0 = (Epsilon_2 + Epsilon_4*su2double(Neighbor_i+1))*StretchingFactor*MeanLambda;
    cte_1 = (Epsilon_2 + Epsilon_4*su2double(Neighbor_j+1))*Area;
//    cte_1 = (Epsilon_2 + Epsilon_4*su2double(Neighbor_j+1))*StretchingFactor*MeanLambda;

    for (iVar = 0; iVar < nVar; iVar++) {
      val_Jacobian_i[iVar][iVar] += fix_factor*cte_0;
      val_Jacobian_j[iVar][iVar] -= fix_factor*cte_1;
    }
  }
}


CSourcePGDMMS::CSourcePGDMMS(unsigned short val_nDim, unsigned short val_nVar, CConfig* config) : CNumerics(val_nDim, val_nVar, config){}


void CSourcePGDMMS::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j,
                                    CConfig* config) {

  const su2double x = Coord_i[0];
  const su2double y = Coord_i[1];


  val_residual[0] = Volume * (exp(x + y)*(20 + cos(x - y) + cos(x + y) + sin(x + y)))/100;

  val_residual[1] = Volume * (exp(x + y)*
                              (200 + 30*cos(y)*
                                         sin(x) +
                               pow(cos(y),2)*
                                   pow(sin(x),2) -
                               10*sin(x)*sin(y) +
                                  cos(x)*
                               (3*pow(cos(y),2)*
                                    sin(x) +
                                   sin(y)*
                                (10 -
                                    sin(x)*sin(y)) +
                                   cos(y)*
                                (30 +
                                 sin(x)*sin(y)))))/
                             100;

  val_residual[2] = Volume * (exp(x + y)*
                              (200 - 10*sin(x)*sin(y) +
                                  pow(cos(x),2)*sin(y)*
                               (3*cos(y) + sin(y)) +
                                  cos(y)*sin(x)*
                               (10 - sin(x)*sin(y)) +
                                  cos(x)*(30*sin(y) +
                                      cos(y)*(30 + sin(x)*sin(y))
                               )))/100;

}

CSourcePTMMS::CSourcePTMMS(unsigned short val_nDim, unsigned short val_nVar, CConfig* config) : CNumerics(val_nDim, val_nVar, config){}


void CSourcePTMMS::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j,
                                    CConfig* config) {

  const su2double k = 1.2*2e-5/18.03e-6;
  const su2double k2 = 4*1000*2e-5*2e-5/ (3*18.03e-6);
  const su2double uf = 10;
  const su2double vf = 10;

  const su2double x = Coord_i[0];
  const su2double y = Coord_i[1];


  val_residual[0] = Volume * (exp(x + y)*(20 + cos(x - y) + cos(x + y) + sin(x + y)))/100;

  val_residual[1] =
      Volume *
      (exp(x + y) *
       (3 * cos(x) * cos(y) * (10 + cos(y) * sin(x)) + pow(10 + cos(y) * sin(x), 2) +
        (10 + cos(y) * sin(x)) * (10 + cos(x) * sin(y)) - sin(x) * sin(y) * (10 + cos(x) * sin(y)) -
        (24 * (-10 + uf - cos(y) * sin(x)) *
         (1 + 0.15 * pow(k * sqrt(pow(10 - uf + cos(y) * sin(x), 2) + pow(10 - vf + cos(x) * sin(y), 2)), 0.687) +
          0.0175 / (1 + 42500. / pow(k * sqrt(pow(10 - uf + cos(y) * sin(x), 2) + pow(10 - vf + cos(x) * sin(y), 2)),
                                     1.16)))) / k2)) / 100;

  val_residual[2] =
      Volume *
      (exp(x + y) *
       (-(sin(x) * (10 + cos(y) * sin(x)) * sin(y)) + 3 * cos(x) * cos(y) * (10 + cos(x) * sin(y)) +
        (10 + cos(y) * sin(x)) * (10 + cos(x) * sin(y)) + pow(10 + cos(x) * sin(y), 2) -
        (24 * (-10 + vf - cos(x) * sin(y)) *
         (1 + 0.15 * pow(k * sqrt(pow(10 - uf + cos(y) * sin(x), 2) + pow(10 - vf + cos(x) * sin(y), 2)), 0.687) +
          0.0175 / (1 + 42500. / pow(k * sqrt(pow(10 - uf + cos(y) * sin(x), 2) + pow(10 - vf + cos(x) * sin(y), 2)),
                                     1.16)))) / k2)) / 100;
}

CSourcePTStagMMS::CSourcePTStagMMS(unsigned short val_nDim, unsigned short val_nVar, CConfig* config) : CNumerics(val_nDim, val_nVar, config){}


void CSourcePTStagMMS::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j,
                                   CConfig* config) {

  const su2double k = 1.2*2e-5/18.03e-6;
  const su2double k2 = 4*1000*2e-5*2e-5/ (3*18.03e-6);
  const su2double uf = 10;
  const su2double vf = 0.0;
  const su2double kk = 7;

  const su2double x = Coord_i[0];
  const su2double y = Coord_i[1];


  val_residual[0] = Volume * (exp(x + y)*kk*(x - y))/100;

  val_residual[1] =
      Volume *
      (exp(x + y) * (pow(kk, 2) * x + pow(kk, 2) * pow(x, 2) - pow(kk, 2) * x * y -
                     (24 * (uf - kk * x) *
                      (1 + 0.15 * pow(k * sqrt(pow(uf - kk * x, 2) + pow(vf + kk * y, 2)), 0.687) +
                       0.0175 / (1 + 42500 / pow(k * sqrt(pow(uf - kk * x, 2) + pow(vf + kk * y, 2)), 1.16)))) /
                         k2)) / 100;

  val_residual[2] =
      Volume *
      (exp(x + y) * (pow(kk, 2) * y - pow(kk, 2) * x * y + pow(kk, 2) * pow(y, 2) -
                     (24 * (vf + kk * y) *
                      (1 + 0.15 * pow(k * sqrt(pow(uf - kk * x, 2) + pow(vf + kk * y, 2)), 0.687) +
                       0.0175 / (1 + 42500 / pow(k * sqrt(pow(uf - kk * x, 2) + pow(vf + kk * y, 2)), 1.16)))) /
                         k2)) / 100;
}