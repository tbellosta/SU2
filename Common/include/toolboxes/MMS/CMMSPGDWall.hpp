
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


#pragma once

#include <cmath>
#include "CVerificationSolution.hpp"

/*!
 * \class CMMSPGDWall
 * \brief Class to define the required data for the manufactured solution of the
          pressureless gas dynamics equations on a quad with a wall BC at the bottom.
 * \author T. Bellosta
 */
class CMMSPGDWall final: public CVerificationSolution {

 public:

  /*!
   * \brief Constructor of the class.
   */
  CMMSPGDWall(void);

  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CMMSPGDWall(unsigned short val_nDim,
                         unsigned short val_nvar,
                         unsigned short val_iMesh,
                         CConfig*       config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMMSPGDWall(void) override;

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const su2double *val_coords,
                   const su2double val_t,
                   su2double       *val_solution) const override;

  void GetPrimitive(const su2double *val_coords,
                            const su2double val_t,
                            su2double       *val_primitive) const override;

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const su2double *val_coords,
                  const su2double val_t,
                  su2double       *val_solution) const override;

  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const su2double *val_coords,
                        const su2double val_t,
                        su2double       *val_source) const override;

  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True, because this is a manufactured solution.
   */
  bool IsManufacturedSolution(void) const override;
};





/*!
 * \class CMMSPGDFar
 * \brief Class to define the required data for the manufactured solution of the
          pressureless gas dynamics equations on a quad with farfield BC.
 * \author T. Bellosta
 */
class CMMSPGDFar final: public CVerificationSolution {

 public:

  /*!
   * \brief Constructor of the class.
   */
  CMMSPGDFar(void);

  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CMMSPGDFar(unsigned short val_nDim,
              unsigned short val_nvar,
              unsigned short val_iMesh,
              CConfig*       config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMMSPGDFar(void) override;

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const su2double *val_coords,
                   const su2double val_t,
                   su2double       *val_solution) const override;

  void GetPrimitive(const su2double *val_coords,
                    const su2double val_t,
                    su2double       *val_primitive) const override;

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const su2double *val_coords,
                  const su2double val_t,
                  su2double       *val_solution) const override;

  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const su2double *val_coords,
                        const su2double val_t,
                        su2double       *val_source) const override;

  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True, because this is a manufactured solution.
   */
  bool IsManufacturedSolution(void) const override;
};














/*!
 * \class CMMSPTFar
 * \brief Class to define the required data for the manufactured solution of the
          particle tracking equations on a quad with farfield BC.
 * \author T. Bellosta
 */
class CMMSPTFar final: public CVerificationSolution {

 public:

  /*!
   * \brief Constructor of the class.
   */
  CMMSPTFar(void);

  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CMMSPTFar(unsigned short val_nDim,
             unsigned short val_nvar,
             unsigned short val_iMesh,
             CConfig*       config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMMSPTFar(void) override;

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const su2double *val_coords,
                   const su2double val_t,
                   su2double       *val_solution) const override;

  void GetPrimitive(const su2double *val_coords,
                    const su2double val_t,
                    su2double       *val_primitive) const override;

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const su2double *val_coords,
                  const su2double val_t,
                  su2double       *val_solution) const override;

  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const su2double *val_coords,
                        const su2double val_t,
                        su2double       *val_source) const override;

  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True, because this is a manufactured solution.
   */
  bool IsManufacturedSolution(void) const override;
};





















/*!
 * \class CMMSPTWall
 * \brief Class to define the required data for the manufactured solution of the
          particle tracking equations on a quad with farfield BC.
 * \author T. Bellosta
 */
class CMMSPTWall final: public CVerificationSolution {

 public:

  /*!
   * \brief Constructor of the class.
   */
  CMMSPTWall(void);

  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CMMSPTWall(unsigned short val_nDim,
            unsigned short val_nvar,
            unsigned short val_iMesh,
            CConfig*       config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMMSPTWall(void) override;

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const su2double *val_coords,
                   const su2double val_t,
                   su2double       *val_solution) const override;

  void GetPrimitive(const su2double *val_coords,
                    const su2double val_t,
                    su2double       *val_primitive) const override;

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const su2double *val_coords,
                  const su2double val_t,
                  su2double       *val_solution) const override;

  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const su2double *val_coords,
                        const su2double val_t,
                        su2double       *val_source) const override;

  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True, because this is a manufactured solution.
   */
  bool IsManufacturedSolution(void) const override;
};
