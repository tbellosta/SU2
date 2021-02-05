/*!
 * \file CMMSPGDWall.hpp
 * \brief Header file for the class CMMSNSUnitQuadSolution.
 *        The implementations are in the <i>CMMSNSUnitQuadSolution.cpp</i> file.
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

#include "../../../include/toolboxes/MMS/CMMSPGDWall.hpp"

CMMSPGDWall::CMMSPGDWall(void) : CVerificationSolution() { }

CMMSPGDWall::CMMSPGDWall(unsigned short val_nDim,
                                               unsigned short val_nVar,
                                               unsigned short val_iMesh,
                                               CConfig*       config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {

  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations on a unit quad. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Solution is being initialized for the " << endl;
    cout << "         manufactured solution  of the" << endl;
    cout << "        particle tracking equations on a unit quad!!!" << endl;
    cout << endl << flush;
  }

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(config->GetTime_Marching() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS PGD Wall case",
                   CURRENT_FUNCTION);

}

CMMSPGDWall::~CMMSPGDWall(void) { }

void CMMSPGDWall::GetBCState(const su2double *val_coords,
                                        const su2double val_t,
                                        su2double       *val_solution) const {

  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CMMSPGDWall::GetPrimitive(const su2double* val_coords, const su2double val_t, su2double* val_primitive) const {

  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  /* Determine the solution for the density, velocity
     components and pressure. */

  const su2double alpha = exp(x + y)/100;

  const su2double u = cos(y)*sin(x) + 10;

  const su2double v = y*(cos(x)*sin(y) + 10);


  /* Compute the conservative variables from the primitive ones.
     Note that the implementation below is valid for both 2D and 3D. */
  val_primitive[0]      = alpha;
  val_primitive[1]      = u;
  val_primitive[2]      = v;

}

void CMMSPGDWall::GetSolution(const su2double *val_coords,
                                         const su2double val_t,
                                         su2double       *val_solution) const {
  GetPrimitive(val_coords, val_t, val_solution);

  for (int iVar = 1; iVar < nVar; ++iVar)
    val_solution[iVar] *= val_solution[0];
}

void CMMSPGDWall::GetMMSSourceTerm(const su2double *val_coords,
                                              const su2double val_t,
                                              su2double       *val_source) const {
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  val_source[0] = (exp(x + y) * (10 * (2 + y) + cos(y) * sin(x) + (1 + y) * cos(x) * (cos(y) + sin(y)))) / 100;

  val_source[1] = (exp(x + y) * (200 + 100 * y + 10 * (3 + y) * cos(y) * sin(x) +
                                 pow(cos(y), 2) * (pow(sin(x), 2) + sin(2 * x)) - 10 * y * sin(x) * sin(y) +
                                 cos(x) * (y * pow(cos(y), 2) * sin(x) + sin(y) * (10 * (1 + y) - y * sin(x) * sin(y)) +
                                           cos(y) * (10 * (2 + y) + (1 + y) * sin(x) * sin(y))))) / 100;

  val_source[2] = (exp(x + y) * y *
       (100 * (3 + y) - 10 * sin(x) * sin(y) + pow(cos(x), 2) * sin(y) * ((1 + 2 * y) * cos(y) + (2 + y) * sin(y)) +
        cos(y) * sin(x) * (10 - sin(x) * sin(y)) +
        cos(x) * (10 * (5 + 2 * y) * sin(y) + cos(y) * (10 + 20 * y + sin(x) * sin(y))))) / 100;
}

bool CMMSPGDWall::IsManufacturedSolution(void) const {
  return true;
}













CMMSPGDFar::CMMSPGDFar(void) : CVerificationSolution() { }

CMMSPGDFar::CMMSPGDFar(unsigned short val_nDim,
                         unsigned short val_nVar,
                         unsigned short val_iMesh,
                         CConfig*       config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {

  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations on a unit quad. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Solution is being initialized for the " << endl;
    cout << "         manufactured solution  of the" << endl;
    cout << "        particle tracking equations on a unit quad!!!" << endl;
    cout << endl << flush;
  }

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(config->GetTime_Marching() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS PGD Far case",
                   CURRENT_FUNCTION);

}

CMMSPGDFar::~CMMSPGDFar(void) { }

void CMMSPGDFar::GetBCState(const su2double *val_coords,
                             const su2double val_t,
                             su2double       *val_solution) const {

  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CMMSPGDFar::GetPrimitive(const su2double* val_coords, const su2double val_t, su2double* val_primitive) const {

  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  /* Determine the solution for the density, velocity
     components and pressure. */

  const su2double alpha = exp(x + y)/100;

  const su2double u = cos(y)*sin(x) + 10;

  const su2double v = cos(x)*sin(y) + 10;


  /* Compute the conservative variables from the primitive ones.
     Note that the implementation below is valid for both 2D and 3D. */
  val_primitive[0]      = alpha;
  val_primitive[1]      = u;
  val_primitive[2]      = v;

}

void CMMSPGDFar::GetSolution(const su2double *val_coords,
                              const su2double val_t,
                              su2double       *val_solution) const {
  GetPrimitive(val_coords, val_t, val_solution);

  for (int iVar = 1; iVar < nVar; ++iVar)
    val_solution[iVar] *= val_solution[0];
}

void CMMSPGDFar::GetMMSSourceTerm(const su2double *val_coords,
                                   const su2double val_t,
                                   su2double       *val_source) const {
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  val_source[0] = (exp(x + y)*(20 + cos(x - y) + cos(x + y) +sin(x + y)))/100;

  val_source[1] = (exp(x + y) *
       (200 + 30 * cos(y) * sin(x) + pow(cos(y), 2) * pow(sin(x), 2) - 10 * sin(x) * sin(y) +
        cos(x) * (3 * pow(cos(y), 2) * sin(x) + sin(y) * (10 - sin(x) * sin(y)) + cos(y) * (30 + sin(x) * sin(y))))) / 100;

  val_source[2] = (exp(x + y)*(200 - 10*sin(x)*sin(y) +pow(cos(x),2)*sin(y)*(3*cos(y) + sin(y)) +
                       cos(y)*sin(x)*(10 - sin(x)*sin(y)) +cos(x)*(30*sin(y) +cos(y)*(30 + sin(x)*sin(y)))))/100;
}

bool CMMSPGDFar::IsManufacturedSolution(void) const {
  return true;
}





























CMMSPTFar::CMMSPTFar(void) : CVerificationSolution() { }

CMMSPTFar::CMMSPTFar(unsigned short val_nDim,
                       unsigned short val_nVar,
                       unsigned short val_iMesh,
                       CConfig*       config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {

  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations on a unit quad. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Solution is being initialized for the " << endl;
    cout << "         manufactured solution  of the" << endl;
    cout << "        particle tracking equations on a unit quad!!!" << endl;
    cout << endl << flush;
  }

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(config->GetTime_Marching() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS PT Far case",
                   CURRENT_FUNCTION);

}

CMMSPTFar::~CMMSPTFar(void) { }

void CMMSPTFar::GetBCState(const su2double *val_coords,
                            const su2double val_t,
                            su2double       *val_solution) const {

  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CMMSPTFar::GetPrimitive(const su2double* val_coords, const su2double val_t, su2double* val_primitive) const {

  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  /* Determine the solution for the density, velocity
     components and pressure. */

  const su2double alpha = exp(x + y)/100;

  const su2double u = cos(y)*sin(x) + 10;

  const su2double v = cos(x)*sin(y) + 10;


  /* Compute the conservative variables from the primitive ones.
     Note that the implementation below is valid for both 2D and 3D. */
  val_primitive[0]      = alpha;
  val_primitive[1]      = u;
  val_primitive[2]      = v;

}

void CMMSPTFar::GetSolution(const su2double *val_coords,
                             const su2double val_t,
                             su2double       *val_solution) const {
  GetPrimitive(val_coords, val_t, val_solution);

  for (int iVar = 1; iVar < nVar; ++iVar)
    val_solution[iVar] *= val_solution[0];
}

void CMMSPTFar::GetMMSSourceTerm(const su2double *val_coords,
                                  const su2double val_t,
                                  su2double       *val_source) const {

  const su2double k = 1.2*2e-5/18.03e-6;
  const su2double k2 = 4*1000*2e-5*2e-5/ (3*18.03e-6);
  const su2double uf = 10;
  const su2double vf = 10;

  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  val_source[0] = (exp(x + y)*(20 + cos(x - y) + cos(x + y) +sin(x + y)))/100;

  val_source[1] = (exp(x + y)*(3*cos(x)*cos(y)*(10 + cos(y)*sin(x)) +pow(10 + cos(y)*sin(x),2) +
                    (10 + cos(y)*sin(x))*(10 + cos(x)*sin(y)) -sin(x)*sin(y)*(10 + cos(x)*sin(y)) -(24*(-10 + uf -
                         cos(y)*sin(x))*(1 +0.15*pow(k*sqrt(pow(10 - uf +cos(y)*sin(x),2) +pow(10 - vf +
                                       cos(x)*sin(y),2)),0.687) +0.0175/(1 +42500./pow(k*sqrt(pow(10 - uf +
                                        cos(y)*sin(x),2) +pow(10 - vf +cos(x)*sin(y),2)),1.16))))/k2))/100;

  val_source[2] = (exp(x + y)*(-(sin(x)*(10 + cos(y)*sin(x))*sin(y)) +3*cos(x)*cos(y)*(10 + cos(x)*sin(y)) +
                    (10 + cos(y)*sin(x))*(10 + cos(x)*sin(y)) +pow(10 + cos(x)*sin(y),2) -(24*(-10 + vf -
                         cos(x)*sin(y))*(1 +0.15*pow(k*sqrt(pow(10 - uf +cos(y)*sin(x),2) +
                                 pow(10 - vf +cos(x)*sin(y),2)),0.687) +0.0175/(1 +42500./pow(k*
                             sqrt(pow(10 - uf +cos(y)*sin(x),2) +pow(10 - vf +cos(x)*sin(y),2)),1.16))))/k2))/100;
}

bool CMMSPTFar::IsManufacturedSolution(void) const {
  return true;
}






























CMMSPTWall::CMMSPTWall(void) : CVerificationSolution() { }

CMMSPTWall::CMMSPTWall(unsigned short val_nDim,
                     unsigned short val_nVar,
                     unsigned short val_iMesh,
                     CConfig*       config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {

  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations on a unit quad. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Solution is being initialized for the " << endl;
    cout << "         manufactured solution  of the" << endl;
    cout << "        particle tracking equations on a unit quad!!!" << endl;
    cout << endl << flush;
  }

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(config->GetTime_Marching() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS PT Wall case",
                   CURRENT_FUNCTION);

}

CMMSPTWall::~CMMSPTWall(void) { }

void CMMSPTWall::GetBCState(const su2double *val_coords,
                           const su2double val_t,
                           su2double       *val_solution) const {

  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CMMSPTWall::GetPrimitive(const su2double* val_coords, const su2double val_t, su2double* val_primitive) const {

  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  /* Determine the solution for the density, velocity
     components and pressure. */

  const su2double alpha = exp(x + y)/100;

  const su2double u = cos(y)*sin(x) + 10;

  const su2double v = y*(cos(x)*sin(y) + 10);


  /* Compute the conservative variables from the primitive ones.
     Note that the implementation below is valid for both 2D and 3D. */
  val_primitive[0]      = alpha;
  val_primitive[1]      = u;
  val_primitive[2]      = v;

}

void CMMSPTWall::GetSolution(const su2double *val_coords,
                            const su2double val_t,
                            su2double       *val_solution) const {
  GetPrimitive(val_coords, val_t, val_solution);

  for (int iVar = 1; iVar < nVar; ++iVar)
    val_solution[iVar] *= val_solution[0];
}

void CMMSPTWall::GetMMSSourceTerm(const su2double *val_coords,
                                 const su2double val_t,
                                 su2double       *val_source) const {

  const su2double k = 1.2*2e-5/18.03e-6;
  const su2double k2 = 4*1000*2e-5*2e-5/ (3*18.03e-6);
  const su2double uf = 10;
  const su2double vf = 10;

  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  val_source[0] = (exp(x + y) * (10 * (2 + y) + cos(y) * sin(x) + (1 + y) * cos(x) * (cos(y) + sin(y)))) / 100;

  val_source[1] =
      (exp(x + y) *
       (2 * cos(x) * cos(y) * (10 + cos(y) * sin(x)) + y * cos(x) * cos(y) * (10 + cos(y) * sin(x)) +
        pow(10 + cos(y) * sin(x), 2) + (10 + cos(y) * sin(x)) * (10 + cos(x) * sin(y)) +
        y * (10 + cos(y) * sin(x)) * (10 + cos(x) * sin(y)) - y * sin(x) * sin(y) * (10 + cos(x) * sin(y)) -
        (24 * (-10 + uf - cos(y) * sin(x)) *
         (1 +
          0.15 * pow(k * sqrt(pow(10 - uf + cos(y) * sin(x), 2) + pow(vf - 10 * y - y * cos(x) * sin(y), 2)),
                       0.687) +
          0.0175 / (1 + 42500.0 / pow(k * sqrt(pow(10 - uf + cos(y) * sin(x), 2) +
                                                 pow(vf - 10 * y - y * cos(x) * sin(y), 2)),
                                        1.16)))) /
            k2)) /
      100;

  val_source[2] =
      (exp(x + y) *
       (-(y * sin(x) * (10 + cos(y) * sin(x)) * sin(y)) + y * cos(x) * cos(y) * (10 + cos(x) * sin(y)) +
        2 * pow(y, 2) * cos(x) * cos(y) * (10 + cos(x) * sin(y)) +
        y * (10 + cos(y) * sin(x)) * (10 + cos(x) * sin(y)) + 2 * y * pow(10 + cos(x) * sin(y), 2) +
        pow(y, 2) * pow(10 + cos(x) * sin(y), 2) -
        (24 * (vf - y * (10 + cos(x) * sin(y))) *
         (1 +
          0.15 * pow(k * sqrt(pow(10 - uf + cos(y) * sin(x), 2) + pow(vf - 10 * y - y * cos(x) * sin(y), 2)),
                       0.687) +
          0.0175 / (1 + 42500.0 / pow(k * sqrt(pow(10 - uf + cos(y) * sin(x), 2) +
                                                 pow(vf - 10 * y - y * cos(x) * sin(y), 2)),
                                        1.16)))) /
            k2)) /
      100;
}

bool CMMSPTWall::IsManufacturedSolution(void) const {
  return true;
}
