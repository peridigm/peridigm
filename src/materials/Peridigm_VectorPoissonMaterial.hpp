//! \file Peridigm_VectorPoissonMaterial.hpp

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

#ifndef PERIDIGM_VECTORPOISSONMATERIAL_HPP
#define PERIDIGM_VECTORPOISSONMATERIAL_HPP

#include "Peridigm_Material.hpp"

namespace PeridigmNS {

  class VectorPoissonMaterial : public Material{
  public:

	//! Constructor.
    VectorPoissonMaterial(const Teuchos::ParameterList & params);

    //! Destructor.
    virtual ~VectorPoissonMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("Vector Poisson"); }

    //! Returns the density of the material.
    virtual double Density() const {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  The Vector Poisson material model is valid only for Peridigm-Albany coupling.\n");
      return 0.0;
    }

    //! Returns the bulk modulus of the material.
    virtual double BulkModulus() const {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  The Vector Poisson material model is valid only for Peridigm-Albany coupling.\n");
      return 0.0;
    }

    //! Returns the shear modulus of the material.
    virtual double ShearModulus() const {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  The Vector Poisson material model is valid only for Peridigm-Albany coupling.\n");
      return 0.0;
    }

    //! Returns a vector of field IDs corresponding to the variables associated with the material.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Returns the requested material property
    //! A dummy method here.
    virtual double lookupMaterialProperty(const std::string keyname) const {return 0.0;}

    //! Initialized data containers and computes weighted volume.
    virtual void
    initialize(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               PeridigmNS::DataManager& dataManager);

    //! Evaluate the internal force.
    virtual void
    computeForce(const double dt,
		 const int numOwnedPoints,
		 const int* ownedIDs,
		 const int* neighborhoodList,
                 PeridigmNS::DataManager& dataManager) const;

  protected:
	
    //! Computes the distance between nodes (a1, a2, a3) and (b1, b2, b3).
    inline double distance(double a1, double a2, double a3,
                           double b1, double b2, double b3) const
      {
        return ( sqrt( (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) + (a3-b3)*(a3-b3) ) );
      }

    // material properties
    double m_horizon;
    double m_coefficient;

    // field spec ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_volumeFieldId;
    int m_modelCoordinatesFieldId;
    int m_coordinatesFieldId;
    int m_forceDensityFieldId;
  };
}

#endif // PERIDIGM_VECTORPOISSONMATERIAL_HPP
