/*!
 * \file CPTVariable.hpp
 * \brief Class for defining the variables of the finite-volume heat equation solver.
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

#include "CVariable.hpp"

/*!
 * \class CPTVariable
 * \brief Class for defining the variables of the finite-volume PT solver.
 * \author O. Burghardt
 * \version 7.0.7 "Blackbird"
 */
class CPTVariable final : public CVariable {
 protected:

  CVectorOfMatrix& Gradient_Reconstruction;  /*!< \brief Reference to the gradient of the primitive variables for MUSCL reconstruction for the convective term */
  CVectorOfMatrix Gradient_Aux;              /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */
  CVectorOfMatrix Gradient_Primitive;              /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */

  MatrixType Primitive;                 /*!< \brief Primitives of the problem. */
  MatrixType Secondary;
  MatrixType Limiter_Primitive;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] volumeFraction - Value of the volume fraction \alpha (initialization value).
   * \param[in] velocity - Value of the particles velocity (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPTVariable(su2double volumeFraction, su2double* velocity, unsigned long npoint, unsigned long ndim,
              unsigned long nvar, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CPTVariable() override = default;

  /*!
   * \brief Get the value of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \param[in] iVar   - Index of the variable.
   * \param[in] iDim   - Index of the dimension.
   * \return Value of the reconstruction variables gradient at a node.
   */
  inline su2double GetGradient_Reconstruction(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const final {
    return Gradient_Reconstruction(iPoint,iVar,iDim);
  }

  /*!
   * \brief Get the value of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \param[in] iVar   - Index of the variable.
   * \param[in] iDim   - Index of the dimension.
   * \param[in] value  - Value of the reconstruction gradient component.
   */
  inline void SetGradient_Reconstruction(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) final {
    Gradient_Reconstruction(iPoint,iVar,iDim) = value;
  }

  /*!
   * \brief Get the array of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \return Array of the reconstruction variables gradient at a node.
   */
  inline su2double **GetGradient_Reconstruction(unsigned long iPoint) final { return Gradient_Reconstruction[iPoint]; }

  /*!
   * \brief Get the reconstruction gradient for primitive variable at all points.
   * \return Reference to variable reconstruction gradient.
   */
  inline CVectorOfMatrix& GetGradient_Reconstruction(void) final { return Gradient_Reconstruction; }


  inline CVectorOfMatrix& GetGradient_Primitive(void) { return Gradient_Primitive; }
  inline const CVectorOfMatrix& GetGradient_Primitive(void) const { return Gradient_Primitive; }

  /*!
   * \brief Get the projected velocity in a unitary vector direction (compressible solver).
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline su2double GetProjVel(unsigned long iPoint, const su2double *val_vector) const final {
    su2double ProjVel = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      ProjVel += Primitive(iPoint,iDim+1)*val_vector[iDim];
    return ProjVel;
  }

  inline su2double GetVelocity(unsigned long iPoint, unsigned long iDim) const final {return Primitive(iPoint,iDim+1); }


  inline void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) final {Primitive(iPoint,iVar) = val_prim;}
  inline void SetPrimitive(unsigned long iPoint, const su2double *val_prim) final {
    for(int iVar = 0; iVar < nVar; iVar++)
      Primitive(iPoint,iVar) = val_prim[iVar];
  }
  inline su2double *GetPrimitive(unsigned long iPoint) final {return Primitive[iPoint]; }

  inline su2double GetPrimitive(unsigned long iPoint, unsigned long iVar) const final {return Primitive(iPoint,iVar); }

  inline const MatrixType& GetPrimitive(void) const { return Primitive; }

  inline su2double *GetLimiter_Primitive(unsigned long iPoint) final { return Limiter_Primitive[iPoint]; }

  inline su2double GetLimiter_Primitive(unsigned long iPoint, unsigned long iVar) const final { return Limiter_Primitive(iPoint,iVar); }

  inline MatrixType& GetLimiter_Primitive(void) {return Limiter_Primitive; }

  inline const MatrixType& GetLimiter_Primitive(void) const {return Limiter_Primitive; }

  inline void SetLimiter_Primitive(unsigned long iPoint, unsigned long iVar, su2double val_value) final {Limiter_Primitive(iPoint,iVar) = val_value;}

  inline su2double **GetGradient_Primitive(unsigned long iPoint) final { return Gradient_Primitive[iPoint]; }

  inline su2double GetGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const final { return Gradient_Primitive(iPoint,iVar,iDim); }

  inline void AddGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) final {
    Gradient_Primitive(iPoint,iVar,iDim) += value;
  }

  inline void SetGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) final {
    Gradient_Primitive(iPoint,iVar,iDim) = value;
  }


  inline su2double GetSecondary(unsigned long iPoint, unsigned long iVar) const final { return Secondary(iPoint,iVar); }

  /*!
   * \brief A virtual member.
   */
  inline void SetSecondary(unsigned long iPoint, unsigned long iVar, su2double val_secondary) final {Secondary(iPoint,iVar) = val_secondary;}

  /*!
   * \brief A virtual member.
   */
  inline void SetSecondary(unsigned long iPoint, const su2double *val_secondary) final {
    for (int iVar = 0; iVar < nSecondaryVar; ++iVar)
      Secondary(iPoint,iVar) = val_secondary[iVar];
  }

  inline su2double *GetSecondary(unsigned long iPoint) final { return Secondary[iPoint]; }

  /*!
 * \brief Set the velocity vector from the old solution.
 * \param[in] val_velocity - Pointer to the velocity.
 */
  inline void SetVelocity_Old(unsigned long iPoint, const su2double *val_velocity) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Solution_Old(iPoint,iDim+1) = val_velocity[iDim]*Solution(iPoint,0);
  }


};
