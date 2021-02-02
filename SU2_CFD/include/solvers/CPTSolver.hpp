/*!
 * \file CPTSolver.hpp
 * \brief Headers of the CPTSolver class
 * \author T.Bellosta
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

#include "CSolver.hpp"
#include "../variables/CPTVariable.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../gradients/computeGradientsGreenGauss.hpp"
#include "../gradients/computeGradientsLeastSquares.hpp"
#include "../limiters/computeLimiters.hpp"

/*!
 * \class CPTSolver
 * \brief Main class for defining the finite-volume PT solver.
 * \author T.Bellosta
 * \version 7.0.7 "Blackbird"
 */
class CPTSolver final : public CSolver {
 protected:
  unsigned short nVarFlow, nMarker, CurrentMesh;
  su2double **HeatFlux, *HeatFlux_per_Marker, *Surface_HF, Total_HeatFlux, AllBound_HeatFlux,
      *AverageT_per_Marker, Total_AverageT, AllBound_AverageT,
      *Primitive, *Primitive_i, *Primitive_j,
      *Surface_Areas, Total_HeatFlux_Areas, Total_HeatFlux_Areas_Monitor;
  su2double ***ConjugateVar, ***InterfaceVar, V_inf;

  su2double **CollectionEfficiency;

  CPTVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

  void clipSolution(void);
  void SetPrimitiveVariables(CGeometry *geometry, CConfig *config);

  void SetCentered_Dissipation_Sensor(CGeometry *geometry, const CConfig *config);

  void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config);

  su2double computeRelaxationTime(CSolver** solver_container, unsigned long iPoint);

  inline su2double uEx(su2double* coord) {
    return 7*coord[0]; //cos(coord[1])*sin(coord[0]) + 10;
  }

  inline su2double vEx(su2double* coord) {
    return -7*coord[1]; //cos(coord[0])*sin(coord[1]) + 10;
  }

  inline su2double aEx(su2double* coord) {
    return exp(coord[0] + coord[1])/100;
  }

 public:

  /*!
   * \brief Constructor of the class.
   */
  CPTSolver(void);

  /*!
   * \brief Constructor of the class.
   */
  CPTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CPTSolver(void) override;

  /*!
   * \brief Restart residual and compute gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry,
                     CSolver **solver_container,
                     CConfig *config,
                     unsigned short iMesh,
                     unsigned short iRKStep,
                     unsigned short RunTime_EqSystem,
                     bool Output) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Postprocessing(CGeometry *geometry,
                      CSolver **solver_container,
                      CConfig *config,
                      unsigned short iMesh) override;

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry,
                   CSolver ***solver,
                   CConfig *config,
                   int val_iter,
                   bool val_update_geo) override;

  /*!
   * \brief Compute the undivided laplacian for the solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute the spatial integration using a centered scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Centered_Residual(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics **numerics_container,
                         CConfig *config,
                         unsigned short iMesh,
                         unsigned short iRKStep) override;
  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) override;

  /*!
   * \brief Compute the viscous residuals for the turbulent equation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Viscous_Residual(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics **numerics_container,
                        CConfig *config,
                        unsigned short iMesh,
                        unsigned short iRKStep) override;

  void Source_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) override;




  /*!
   * \brief Update the solution using an implicit solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitEuler_Iteration(CGeometry *geometry,
                               CSolver **solver_container,
                               CConfig *config) override;

  /*!
   * \brief Update the solution using an explicit solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ExplicitEuler_Iteration(CGeometry *geometry,
                               CSolver **solver_container,
                               CConfig *config) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry,
                    CSolver **solver_container,
                    CConfig *config,
                    unsigned short iMesh,
                    unsigned long Iteration) override;

  /*!
   * \brief Set the initial condition for the FEM structural problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry,
                           CSolver ***solver_container,
                           CConfig *config,
                           unsigned long TimeIter) override;

  /*!
   * \brief Set the total residual adding the term that comes from the Dual Time-Stepping Strategy.
   * \param[in] geometry - Geometric definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void SetResidual_DualTime(CGeometry *geometry,
                            CSolver **solver_container,
                            CConfig *config,
                            unsigned short iRKStep,
                            unsigned short iMesh,
                            unsigned short RunTime_EqSystem) override;


  void BC_Far_Field(CGeometry *geometry,
                    CSolver **solver_container,
                    CNumerics *conv_numerics,
                    CNumerics *visc_numerics,
                    CConfig *config,
                    unsigned short val_marker) final;

  void BC_HeatFlux_Wall(CGeometry *geometry,
                    CSolver **solver_container,
                    CNumerics *conv_numerics,
                    CNumerics *visc_numerics,
                    CConfig *config,
                    unsigned short val_marker) final;

  void BC_Euler_Wall(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) final;

  inline su2double GetCollectionEfficiency(unsigned short val_marker, unsigned long val_vertex) const final {
    return CollectionEfficiency[val_marker][val_vertex];
  }

  void computeCollectionEfficiency(CGeometry *geometry,
                                   CSolver **solver_container,
                                   CConfig *config,
                                   unsigned short iMesh);


  inline void SetPrimitive_Gradient_GG(CGeometry *geometry,
                                               const CConfig *config,
                                               bool reconstruction = false) final {

    const auto& primitives = nodes->GetPrimitive();
    auto& gradient = reconstruction ? nodes->GetGradient_Reconstruction() : nodes->GetGradient_Primitive();

    computeGradientsGreenGauss(this, PRIMITIVE_GRADIENT, PERIODIC_PRIM_GG, *geometry, *config, primitives, 0,
                               nPrimVar, gradient);

  }

  inline void SetPrimitive_Gradient_LS(CGeometry *geometry,
                                       const CConfig *config,
                                       bool reconstruction = false) final {

    /*--- Set a flag for unweighted or weighted least-squares. ---*/
    bool weighted;

    if (reconstruction)
      weighted = (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES);
    else
      weighted = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);

    const auto& primitives = nodes->GetPrimitive();
    auto& rmatrix = nodes->GetRmatrix();
    auto& gradient = reconstruction ? nodes->GetGradient_Reconstruction() : nodes->GetGradient_Primitive();
    PERIODIC_QUANTITIES kindPeriodicComm = weighted ? PERIODIC_PRIM_LS : PERIODIC_PRIM_ULS;

    computeGradientsLeastSquares(this, PRIMITIVE_GRADIENT, kindPeriodicComm, *geometry, *config, weighted, primitives, 0,
                                 nPrimVar, gradient, rmatrix);

  }

  inline void SetPrimitive_Limiter(CGeometry *geometry, const CConfig *config) final {

    auto kindLimiter = static_cast<ENUM_LIMITER>(config->GetKind_SlopeLimit_PT());
    const auto& primitives = nodes->GetPrimitive();
    const auto& gradient = nodes->GetGradient_Reconstruction();
    auto& primMin = nodes->GetSolution_Min();
    auto& primMax = nodes->GetSolution_Max();
    auto& limiter = nodes->GetLimiter_Primitive();

    computeLimiters(kindLimiter, this, PRIMITIVE_LIMITER, PERIODIC_LIM_PRIM_1, PERIODIC_LIM_PRIM_2, *geometry, *config, 0,
                    nPrimVar, primitives, gradient, primMin, primMax, limiter);

  }


};