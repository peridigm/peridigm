/*! \file Peridigm_Material.hpp */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#ifndef PERIDIGM_MATERIAL_HPP
#define PERIDIGM_MATERIAL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <vector>
#include <string>
#include <float.h>
#include "Peridigm_DataManager.hpp"
#include "Peridigm_SerialMatrix.hpp"
#include "Peridigm_ScratchMatrix.hpp"

namespace PeridigmNS {

  //! Base class defining the Peridigm material model interface.
  class Material{

  public:

    //! Standard constructor.
    Material(const Teuchos::ParameterList & params) : m_finiteDifferenceProbeLength(DBL_MAX) {
      if(params.isParameter("Finite Difference Probe Length"))
      m_finiteDifferenceProbeLength = params.get<double>("Finite Difference Probe Length");
    }

    //! Destructor.
    virtual ~Material(){}

    //! Return name of material type
    virtual std::string Name() const = 0;

    //! Returns the density of the material.
    virtual double Density() const = 0;

    //! Returns the bulk modulus of the material.
    virtual double BulkModulus() const = 0;

    //! Returns the shear modulus of the material.
    virtual double ShearModulus() const = 0;

	//! Returns material property value for a given key
	// Only implemented for multiphysics elastic material
	virtual double lookupMaterialProperty(const std::string keyname) const {
      std::string errorMsg = "**Error, Material::lookupMaterialProperty() called for ";
      errorMsg += Name();
      errorMsg += " but this function is not implemented.\n";
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, errorMsg);
      return 0.0;
    }

    //! Returns a vector of field IDs corresponding to the variables associated with the material.
    virtual std::vector<int> FieldIds() const = 0;

    //! Returns a vector of field IDs that need to be synchronized across block boundaries and MPI boundaries after initialize().
    virtual std::vector<int> FieldIdsForSynchronizationAfterInitialize() const {
      std::vector<int> empty;
      return empty;
    }

    //! Returns a vector of field IDs that need to be synchronized across block boundaries and MPI boundaries after precompute().
    virtual std::vector<int> FieldIdsForSynchronizationAfterPrecompute() const {
      std::vector<int> empty;
      return empty;
    }

    //! Initialize the material model.
    virtual void
    initialize(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               PeridigmNS::DataManager& dataManager) {}

    //! Run calculations at each time step prior to evaulating the internal force.
    virtual void
    precompute(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               PeridigmNS::DataManager& dataManager) const {}

    //! Evaluate the internal force.
    virtual void
    computeForce(const double dt,
                 const int numOwnedPoints,
                 const int* ownedIDs,
                 const int* neighborhoodList,
                 PeridigmNS::DataManager& dataManager) const {};

    //! Compute the divergence of the flux (for diffusion models).
    virtual void
    computeFluxDivergence(const double dt,
                          const int numOwnedPoints,
                          const int* ownedIDs,
                          const int* neighborhoodList,
                          PeridigmNS::DataManager& dataManager) const {}

    /// \enum JacobianType
    /// \brief Whether to compute the full tangent stiffness matrix or just its block diagonal entries
    ///
    /// The Peridigm Material base class provides a computeJacobian method that all materials inherit
    /// to compute the tangent stiffness matrix. This base class uses a finite difference method (it provides both
    /// forward and centered) to numerically approximate the jacobian. Derived classes may override this method
    /// to compute the jacobian via another approach (for example, automatic differentiation) or may simply
    /// inherit and use the finite difference Jacobian, which will work for all derived mateiral classes.
    ///
    /// The default behavior of this
    /// method is to compute the full tangent stiffness matrix. However, is it occasionally useful to compute
    /// and store only the block diagonal entries -- specifically, the only the entries of the matrix that
    /// describe interactions between different dofs for an individual node. This manifests as a block-diagonal
    /// with each block of dimension 3x3, and one block for each node. This block-diagonal matrix is useful,
    /// for example, as a preconditioner for the full Jacobian matrix. The 3x3 blocks are also the "P" matrices
    /// used by the ComputeStabilityIndex compute class. See that class for more information.
    ///
    /// \note The default behavior is to compute the full tangent stiffness matrix. This enum is useful to only
    /// if you need to efficiently compute only the block diagonal entries of the full tangent stiffness matrix.
    enum JacobianType { UNDEFINED=0, NONE=1, FULL_MATRIX=2, BLOCK_DIAGONAL=3 };

    //! Evaluate the jacobian.
    virtual void
    computeJacobian(const double dt,
                    const int numOwnedPoints,
                    const int* ownedIDs,
                    const int* neighborhoodList,
                    PeridigmNS::DataManager& dataManager,
                    PeridigmNS::SerialMatrix& jacobian,
                    PeridigmNS::Material::JacobianType jacobianType = PeridigmNS::Material::FULL_MATRIX) const;

    //! Compute stored elastic energy density
    virtual void
    computeStoredElasticEnergyDensity(const double dt,
                                      const int numOwnedPoints,
                                      const int* ownedIDs,
                                      const int* neighborhoodList,
                                      PeridigmNS::DataManager& dataManager) const {}

    //! Compute the bulk modulus given any two elastic constants from among:  bulk modulus, shear modulus, Young's modulus, Poisson's ratio.
    double calculateBulkModulus(const Teuchos::ParameterList & params) const;

    //! Compute the shear modulus given any two elastic constants from among:  bulk modulus, shear modulus, Young's modulus, Poisson's ratio.
    double calculateShearModulus(const Teuchos::ParameterList & params) const;

    enum FiniteDifferenceScheme { FORWARD_DIFFERENCE=0, CENTRAL_DIFFERENCE=1 };

  protected:

    //! Evaluate the jacobian via finite difference (probing)
    void
    computeFiniteDifferenceJacobian(const double dt,
                                    const int numOwnedPoints,
                                    const int* ownedIDs,
                                    const int* neighborhoodList,
                                    PeridigmNS::DataManager& dataManager,
                                    PeridigmNS::SerialMatrix& jacobian,
                                    FiniteDifferenceScheme finiteDifferenceScheme,
                                    PeridigmNS::Material::JacobianType jacobianType = PeridigmNS::Material::FULL_MATRIX) const;

    //! Scratch matrix.
    mutable ScratchMatrix scratchMatrix;

    //! Finite-difference probe length
    double m_finiteDifferenceProbeLength;

  private:

    //! Default constructor with no arguments, private to prevent use.
    Material(){}
  };
}

#endif // PERIDIGM_MATERIAL_HPP
