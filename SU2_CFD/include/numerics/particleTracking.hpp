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
      **Proj_Jac_i, **Proj_Jac_j, *Conservatives_i, *Conservatives_j, *IntermediateState;
  bool implicit, dynamic_grid;
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

   void GetProjFluxPT(const su2double* VolFraction, const su2double* Vel, const su2double* Pi, su2double* ProjFlux,
                      const su2double* Norm) const;

   void GetProjFluxJacobianPT(const su2double* VolFraction, const su2double* Vel, const su2double* Pi,
                              const su2double* Norm, su2double** ProjJac) const;
                              
  
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


/*!
 * \class CUpwGodunov_PT
 * \brief Class for doing a scalar upwind solver for the heat convection equation.
 * \ingroup ConvDiscr
 * \author T. Bellosta
 * \version 7.0.7 "Blackbird"
 */
class CUpwGodunov_PT : public CConv_PT {

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwGodunov_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  ~CUpwGodunov_PT() override = default;


  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) final;
};


/*!
 * \class CUpwStegWarm_PT
 * \brief Class for doing a scalar upwind solver for the heat convection equation.
 * \ingroup ConvDiscr
 * \author T. Bellosta
 * \version 7.0.7 "Blackbird"
 */
class CUpwStegWarm_PT : public CConv_PT {

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwStegWarm_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  ~CUpwStegWarm_PT() override = default;


  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) final;
};




/*!
 * \class CUpwGodunov_PT
 * \brief Class for doing a scalar upwind solver for the heat convection equation.
 * \ingroup ConvDiscr
 * \author T. Bellosta
 * \version 7.0.7 "Blackbird"
 */
class CUpwFDS_PT : public CConv_PT {

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwFDS_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  ~CUpwFDS_PT() override = default;


  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) final;
};


/*!
 * \class CUpwHLLC_PT
 * \brief Class for doing a scalar upwind solver for the heat convection equation.
 * \ingroup ConvDiscr
 * \author T. Bellosta
 * \version 7.0.7 "Blackbird"
 */
class CUpwHLLC_PT : public CConv_PT {
 protected:

  su2double a2;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwHLLC_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  ~CUpwHLLC_PT() override = default;


  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) final;
};


class CSourceDrag : public CNumerics {
 private:
  su2double *flowPrimitive;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceDrag(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  ~CSourceDrag() override = default;


  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) final;
};


/*!
 * \class CAvgGradCorrected_PT
 * \brief Class for computing viscous term using average of gradients with correction (heat equation).
 * \ingroup ViscDiscr
 * \author O. Burghardt.
 * \version 7.0.7 "Blackbird"
 */
class CAvgGradCorrected_PT : public CNumerics {
 protected:
  su2double ** Mean_GradPrimVar;
  su2double *Proj_Mean_GradPrimVar_Kappa, *Proj_Mean_GradPrimVar_Edge, *Proj_Mean_GradPrimVar_Corrected;
  su2double *Edge_Vector;
  bool implicit;
  su2double dist_ij_2, proj_vector_ij, viscosity_Mean;
  unsigned short iVar, iDim;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_PT(void) override;

  /*!
   * \brief Compute the viscous heat residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) override;
};


/*!
 * \class CAvgGrad_PT
 * \brief Class for computing viscous term using average of gradients with correction (heat equation).
 * \ingroup ViscDiscr
 * \author O. Burghardt.
 * \version 7.0.7 "Blackbird"
 */
class CAvgGrad_PT : public CAvgGradCorrected_PT {
 protected:


 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_PT(void) override = default;

  /*!
   * \brief Compute the viscous heat residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) override;
};


class CCent_PT : public CConv_PT {
 private:
  unsigned short iVar, jVar, iDim;
  su2double *MeanVelocity, MeanDensity, *Diff_U, *Diff_Lapl;
  su2double MeanLambda, StrechingFactor, Param_p, StretchingFactor,
      sc2, sc4, Epsilon_2, Epsilon_4, Param_Kappa_2, Param_Kappa_4, cte_0, cte_1, fix_factor;

  void DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j);

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCent_PT(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  ~CCent_PT() override;


  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) final;

};
