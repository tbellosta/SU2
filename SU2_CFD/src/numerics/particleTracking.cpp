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
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

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

  IntermediateState = nullptr;

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

  delete [] IntermediateState;
}

void CConv_PT::GetProjFluxPT(const su2double* VolFraction, const su2double* Vel, const su2double* Pi,
                             su2double* ProjFlux, const su2double* Norm) const {

  su2double u = Vel[0];
  su2double v = Vel[1];
  su2double P = (*Pi);
  su2double q = GeometryToolbox::DotProduct(nDim,Vel,Norm);
  su2double aq = (*VolFraction)*q;


  ProjFlux[0] = aq;
  ProjFlux[1] = aq*u + P*Norm[0];
  ProjFlux[2] = aq*v + P*Norm[1];
  if (nDim == 3) ProjFlux[3] = aq*Vel[2]+P*Norm[2];

}

void CConv_PT::GetProjFluxJacobianPT(const su2double* VolFraction, const su2double* Vel, const su2double* Pi,
                                     const su2double* Norm, su2double** ProjJac) const {

  su2double alpha = *VolFraction;
  su2double u = Vel[0];
  su2double v = Vel[1];
  su2double P = (*Pi);
  su2double q = GeometryToolbox::DotProduct(nDim,Vel,Norm);
  su2double aq = (*VolFraction)*q;

  const su2double a = RelaxationFactor;

  ProjJac[0][0] = 0;                     ProjJac[0][1] = Norm[0];                     ProjJac[0][2] = Norm[1];
  ProjJac[1][0] = -u*q-P*Norm[0]/alpha;  ProjJac[1][1] = q+u*Norm[0];                 ProjJac[1][2] = u*Norm[1];
  ProjJac[2][0] = -v*q-P*Norm[1]/alpha;  ProjJac[2][1] = v*Norm[0];                   ProjJac[2][2] = q+v*Norm[1];

  if (nDim == 3) SU2_MPI::Error("Jacobian not implemented for 3D cases.", CURRENT_FUNCTION);
}


CUpwRusanov_PT::CUpwRusanov_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
    CConv_PT(val_nDim, val_nVar, config) { }


void CUpwRusanov_PT::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                   su2double **val_Jacobian_j, CConfig *config) {
#define SMALL 1e-50

  const su2double gamma = 1.4;

  su2double projVel_i = 0;
  su2double projVel_j = 0;
  su2double lambda, E_i, E_j, MagVel2_i, MagVel2_j, c_i, c_j;

  Density_i = V_i[0];
  Density_j = V_j[0];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    projVel_i += Velocity_i[iDim]*Normal[iDim];
    projVel_j += Velocity_j[iDim]*Normal[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Pressure_j = V_j[nDim+1];

  MagVel2_i = GeometryToolbox::SquaredNorm(nDim,Velocity_i);
  MagVel2_j = GeometryToolbox::SquaredNorm(nDim,Velocity_j);

  E_i = 0.5*MagVel2_i + Pressure_i / (Density_i*(gamma-1));
  E_j = 0.5*MagVel2_j + Pressure_j / (Density_j*(gamma-1));

  c_i = (Pressure_i < PT_EPS) ? PT_EPS : sqrt(gamma*Pressure_i/Density_i);
  c_j = (Pressure_j < PT_EPS) ? PT_EPS : sqrt(gamma*Pressure_j/Density_j);

  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;

  for (iDim = 0; iDim < nDim; ++iDim) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  Conservatives_i[nVar-1] = Density_i*E_i;
  Conservatives_j[nVar-1] = Density_j*E_j;

  lambda = fmax(fabs(projVel_i)+c_i,fabs(projVel_j)+c_j);

  GetProjFluxPT(&Density_i, Velocity_i, &E_i, Proj_Flux_i, Normal);
  GetProjFluxPT(&Density_j, Velocity_j, &E_j, Proj_Flux_j, Normal);

  for (int iVar = 0; iVar < nVar; ++iVar) {
    val_residual[iVar] = 0.5*(Proj_Flux_i[iVar]+Proj_Flux_j[iVar]) - 0.5*lambda*(Conservatives_j[iVar] - Conservatives_i[iVar]);
  }


  if (implicit) {
    GetProjFluxJacobianPT(&Density_i, Velocity_i, &E_i, Normal, Proj_Jac_i);
    GetProjFluxJacobianPT(&Density_j, Velocity_j, &E_j, Normal, Proj_Jac_j);

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
    CConv_PT(val_nDim, val_nVar, config) {

  IntermediateState = new su2double [nVar];

}


void CUpwGodunov_PT::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                     su2double **val_Jacobian_j, CConfig *config) {
#define MAXNDIM 3

  const su2double a = RelaxationFactor;

  su2double projVel_i = 0;
  su2double projVel_j = 0;
  su2double Area, sL, sM, sR, pStar, aStar;
  su2double UnitNormal[MAXNDIM];

  Area = GeometryToolbox::Norm(nDim,Normal);

  Density_i = V_i[0];
  Density_j = V_j[0];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    UnitNormal[iDim] = Normal[iDim] / Area;
    projVel_i += Velocity_i[iDim]*UnitNormal[iDim];
    projVel_j += Velocity_j[iDim]*UnitNormal[iDim];
  }

  Pressure_i = 0.0;
  Pressure_j = 0.0;

  sL = projVel_i - (a/Density_i);
  sM = ( projVel_i + projVel_j) / 2;
  sR = projVel_j + (a/Density_j);

  pStar = a*(projVel_i-projVel_j) / 2;

  if (sL >= 0) {
    GetProjFluxPT(&Density_i, Velocity_i, &Pressure_i, val_residual, UnitNormal);
  } else if (sM > 0) {

    aStar = a*Density_i / (Density_i*(sM-projVel_i) + a);

    IntermediateState[0] = aStar;
    for (int iDim = 0; iDim < nDim; ++iDim) IntermediateState[iDim+1] = Velocity_i[iDim] + (sM-projVel_i)*UnitNormal[iDim];

    GetProjFluxPT(&aStar, &(IntermediateState[1]), &pStar, val_residual, UnitNormal);

  } else if (sR > 0) {

    aStar = a*Density_j / (Density_j*(projVel_j-sM) + a);

    IntermediateState[0] = aStar;
    for (int iDim = 0; iDim < nDim; ++iDim) IntermediateState[iDim+1] = Velocity_j[iDim] + (sM-projVel_j)*UnitNormal[iDim];

    GetProjFluxPT(&aStar, &(IntermediateState[1]), &pStar, val_residual, UnitNormal);

  } else {
    GetProjFluxPT(&Density_j, Velocity_j, &Pressure_j, val_residual, UnitNormal);
  }

  for (int iVar = 0; iVar < nVar; ++iVar) val_residual[iVar] *= Area;

  if (implicit) {

    for (int iVar = 0; iVar < nVar; ++iVar) {
      for (int jVar = 0; jVar < nVar; ++jVar) {
        val_Jacobian_i[iVar][jVar] = 0.0;
        val_Jacobian_j[iVar][jVar] = 0.0;
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


  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    projVel_i += Velocity_i[iDim]*Normal[iDim];
    projVel_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;

  for (iDim = 0; iDim < nDim; ++iDim) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }

  uRoe = (sqrt(Density_i)*projVel_i + sqrt(Density_j)*projVel_j) / (sqrt(Density_i) + sqrt(Density_j));
  uRoe /= Area;
  uRoe = fabs(uRoe);

  su2double eps = 1e-2;

  if (uRoe < eps) uRoe = 0.5*(uRoe*uRoe/eps + eps);

  GetProjFluxPT(&Density_i, Velocity_i, nullptr, Proj_Flux_i, Normal);
  GetProjFluxPT(&Density_j, Velocity_j, nullptr, Proj_Flux_j, Normal);

  for (int iVar = 0; iVar < nVar; ++iVar) {
    val_residual[iVar] = 0.5 * ((Proj_Flux_i[iVar] + Proj_Flux_j[iVar]) - uRoe * (Conservatives_j[iVar] - Conservatives_i[iVar]) * Area);
  }

  if (implicit) {
    GetProjFluxJacobianPT(&Density_i, Velocity_i, nullptr, Normal, Proj_Jac_i);
    GetProjFluxJacobianPT(&Density_j, Velocity_j, nullptr, Normal, Proj_Jac_j);

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




CUpwHLLC_PT::CUpwHLLC_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
    CConv_PT(val_nDim, val_nVar, config) {

  IntermediateState = new su2double[nVar];
  a2 = config->GetParticle_Size()*STANDARD_GRAVITY;

}

void CUpwHLLC_PT::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  su2double projVel_i = 0;
  su2double projVel_j = 0;
  su2double uStar, aStar, sL, sM, sR, qL, qR, ratio, aSplit, diffSpeed;

//  const su2double a = sqrt(a2);

  Density_i = V_i[0];
  Density_j = V_j[0];

//  const su2double aL = sqrt(STANDARD_GRAVITY*Density_i);
//  const su2double aR = sqrt(STANDARD_GRAVITY*Density_j);

  const su2double aL = 100;
  const su2double aR = 100;

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;

  for (iDim = 0; iDim < nDim; ++iDim) {
    UnitNormal[iDim] = Normal[iDim] / Area;
    projVel_i += Velocity_i[iDim]*UnitNormal[iDim];
    projVel_j += Velocity_j[iDim]*UnitNormal[iDim];

    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }

  for (int iVar = 0; iVar < nVar; ++iVar) val_residual[iVar] = 0.0;

  aStar =
      0.5 * (Density_i + Density_j) - 0.25 * (projVel_j - projVel_i) * (Density_i + Density_j) / (aL + aR);
  uStar =
      0.5 * (projVel_i + projVel_j) - 0.25 * (Density_j - Density_i) * (aL + aR) / (Density_i + Density_j);

//  aStar = sqrt(Density_i*Density_j)*exp((projVel_i-projVel_j)/(2*a));
//  uStar = 0.5*(projVel_i+projVel_j) + 0.5*a*log(Density_i/Density_j);

//  qL = (aStar > Density_i) ? sqrt(0.5*(aStar+Density_i)*aStar/pow(Density_i,2)) : 1.0;
//  qR = (aStar > Density_j) ? sqrt(0.5*(aStar+Density_j)*aStar/pow(Density_j,2)) : 1.0;

  //Gave errors for some unknown reason (?)
  if(aStar > Density_i){
    qL = sqrt(aStar/Density_i);
  }else{
    qL = 1.0;
  }
  if(aStar > Density_j){
    qR = sqrt(aStar/Density_j);
  }else{
    qR = 1.0;
  }

  sL = projVel_i - aL * qL;
  sR = projVel_j + aR * qR;
  sM = (sL * Density_j * (projVel_j - sR) - sR * Density_i * (projVel_i - sL));
  if (sM!=0.0) sM /= (Density_j * (projVel_j - sR) - Density_i * (projVel_i - sL));
//  sM = uStar;

  if (sL >= 0) {
    val_residual[0] = Density_i * projVel_i;
    for (int iVar = 1; iVar < nVar; ++iVar)
//      val_residual[iVar] = Density_i * projVel_i * Velocity_i[iVar - 1] +  0.5 * STANDARD_GRAVITY * pow(Density_i,2) * UnitNormal[iVar - 1];
      val_residual[iVar] = Density_i * projVel_i * Velocity_i[iVar - 1] +  Density_i * a2 * UnitNormal[iVar - 1];

    aSplit = Density_i;

  } else if (sM >= 0) {

    diffSpeed = sL- sM;
    ratio = (sL - projVel_i) / diffSpeed;

    IntermediateState[0] = ratio * Density_i;
    for (iDim = 0; iDim < nDim; ++iDim) {
//      IntermediateState[iDim+1] = (Density_i*Velocity_i[iDim]*(sL-projVel_i) + 0.5 * STANDARD_GRAVITY * (pow(IntermediateState[0],2) - pow(Density_i,2))*UnitNormal[iDim]) / diffSpeed;
      IntermediateState[iDim+1] = (Density_i*Velocity_i[iDim]*(sL-projVel_i) + a2*(IntermediateState[0]-Density_i)*UnitNormal[iDim]) / diffSpeed;
    }

    val_residual[0] = Density_i * projVel_i;
    for (int iVar = 1; iVar < nVar; ++iVar)
//      val_residual[iVar] = Density_i * projVel_i * Velocity_i[iVar - 1] +  0.5 * STANDARD_GRAVITY * pow(Density_i,2) * UnitNormal[iVar - 1];
      val_residual[iVar] = Density_i * projVel_i * Velocity_i[iVar - 1] +  Density_i * a2 * UnitNormal[iVar - 1];

    val_residual[0] += sL * (IntermediateState[0] - Density_i);
    for (int iVar = 1; iVar < nVar; ++iVar) {
      val_residual[iVar] += sL * (IntermediateState[iVar] - Conservatives_i[iVar]);
    }

    aSplit = IntermediateState[0];

  } else if (sR >= 0) {

    diffSpeed = sR- sM;
    ratio = (sR - projVel_j) / diffSpeed;

    IntermediateState[0] = ratio * Density_j;
    for (iDim = 0; iDim < nDim; ++iDim) {
//      IntermediateState[iDim+1] = (Density_j*Velocity_j[iDim]*(sR-projVel_j) + 0.5 * STANDARD_GRAVITY * (pow(IntermediateState[0],2) - pow(Density_j,2))*UnitNormal[iDim]) / diffSpeed;
      IntermediateState[iDim+1] = (Density_j*Velocity_j[iDim]*(sR-projVel_j) + a2*(IntermediateState[0]-Density_j)*UnitNormal[iDim]) / diffSpeed;
    }

    val_residual[0] = Density_j * projVel_j;
    for (int iVar = 1; iVar < nVar; ++iVar)
//      val_residual[iVar] = Density_j * projVel_j * Velocity_j[iVar - 1] +  0.5 * STANDARD_GRAVITY * pow(Density_j,2) * UnitNormal[iVar - 1];
      val_residual[iVar] = Density_j * projVel_j * Velocity_j[iVar - 1] +  Density_j * a2 * UnitNormal[iVar - 1];

    val_residual[0] += sR * (IntermediateState[0] - Density_j);
    for (int iVar = 1; iVar < nVar; ++iVar) {
      val_residual[iVar] += sR * (IntermediateState[iVar] - Conservatives_j[iVar]);
    }

    aSplit = IntermediateState[0];

  } else {
    val_residual[0] = Density_j * projVel_j;
    for (int iVar = 1; iVar < nVar; ++iVar)
//      val_residual[iVar] = Density_j * projVel_j * Velocity_j[iVar - 1] +  0.5 * STANDARD_GRAVITY * pow(Density_j,2) * UnitNormal[iVar - 1];
      val_residual[iVar] = Density_j * projVel_j * Velocity_j[iVar - 1] +  Density_j * a2 * UnitNormal[iVar - 1];

    aSplit = Density_j;
  }

  if (aSplit <= 0) {throw std::runtime_error("error");}

  for (int iVar = 0; iVar < nVar; ++iVar) val_residual[iVar] *= Area;

//  for (iDim = 0; iDim < nDim; ++iDim) val_residual[iDim+1] -= a2*aSplit*Normal[iDim];

}


void CSourceDrag::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j,
                                  CConfig* config) {

  su2double *uFlow = &V_i[1];

//  su2double Uf[2] = {4.0, 1.0};
//  uFlow = Uf;

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

  const su2double gamma = 1.4;

  su2double projVel_i = 0;
  su2double projVel_j = 0;

  su2double MagVel2_i, MagVel2_j, E_i, E_j;

  Density_i = V_i[0];
  Density_j = V_j[0];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    projVel_i += Velocity_i[iDim]*Normal[iDim];
    projVel_j += Velocity_j[iDim]*Normal[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Pressure_j = V_j[nDim+1];

  MagVel2_i = GeometryToolbox::SquaredNorm(nDim,Velocity_i);
  MagVel2_j = GeometryToolbox::SquaredNorm(nDim,Velocity_j);

  E_i = 0.5*MagVel2_i + Pressure_i / (Density_i*(gamma-1));
  E_j = 0.5*MagVel2_j + Pressure_j / (Density_j*(gamma-1));

  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;

  for (iDim = 0; iDim < nDim; ++iDim) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  Conservatives_i[nVar-1] = Density_i*E_i;
  Conservatives_j[nVar-1] = Density_j*E_j;

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

  GetProjFluxPT(&Density_i, Velocity_i, nullptr, Proj_Flux_i, Normal);
  GetProjFluxPT(&Density_j, Velocity_j, nullptr, Proj_Flux_j, Normal);

  /*--- Jacobians of the inviscid flux, scale = 0.5 because ProjFlux ~ 0.5*(fc_i+fc_j)*Normal ---*/

  if (implicit) {
    GetProjFluxJacobianPT(&Density_i, Velocity_i, nullptr, Normal, Proj_Jac_i);
    GetProjFluxJacobianPT(&Density_j, Velocity_j, nullptr, Normal, Proj_Jac_j);
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
