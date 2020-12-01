/*!
 * \file heat.hpp
 * \brief Delarations of numerics classes for particle tracking problems.
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

#pragma once

#include "CNumerics.hpp"

class CConv_PT : public CNumerics {
 protected:

  su2double *Velocity_i, *Velocity_j, *Proj_Flux_i, *Proj_Flux_j,
      **Proj_Jac_i, **Proj_Jac_j;
  bool implicit, dynamic_grid;
//  su2double q_ij, a0, a1;
  unsigned short iDim;

 public:

  /*!
 * \brief Constructor of the class.
 * \param[in] val_nDim - Number of dimensions of the problem.
 * \param[in] val_nVar - Number of variables of the problem.
 * \param[in] config - Definition of the particular problem.
 */
  CConv_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CConv_PT() override;

  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
   void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) override = 0;

   void GetProjFluxPT(const su2double * VolFraction, const su2double * Vel, const su2double *Normal,
                      su2double* ProjFlux) const;

  void GetProjFluxJacobianPT(const su2double * VolFraction, const su2double * Vel, const su2double *Normal,
                              su2double** ProjJac) const;

};


/*!
 * \class CUpwRusanov_PT
 * \brief Class for doing a scalar upwind solver for the heat convection equation.
 * \ingroup ConvDiscr
 * \author T. Bellosta
 * \version 7.0.7 "Blackbird"
 */
class CUpwRusanov_PT : public CConv_PT {

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRusanov_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  ~CUpwRusanov_PT() override = default;


  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) final;
};
