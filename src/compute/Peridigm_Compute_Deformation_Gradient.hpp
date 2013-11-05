/*! \file Peridigm_Compute_Deformation_Gradient.hpp */

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

#ifdef COMPUTE_CLASS

ComputeClass(Deformation_GradientXX,Compute_Deformation_Gradient)
ComputeClass(Deformation_GradientXY,Compute_Deformation_Gradient)
ComputeClass(Deformation_GradientXZ,Compute_Deformation_Gradient)
ComputeClass(Deformation_GradientYX,Compute_Deformation_Gradient)
ComputeClass(Deformation_GradientYY,Compute_Deformation_Gradient)
ComputeClass(Deformation_GradientYZ,Compute_Deformation_Gradient)
ComputeClass(Deformation_GradientZX,Compute_Deformation_Gradient)
ComputeClass(Deformation_GradientZY,Compute_Deformation_Gradient)
ComputeClass(Deformation_GradientZZ,Compute_Deformation_Gradient)

#else

#ifndef PERIDIGM_COMPUTE_DEFORMATION_GRADIENT_HPP
#define PERIDIGM_COMPUTE_DEFORMATION_GRADIENT_HPP

#include "Peridigm_Compute.hpp"

namespace PeridigmNS {

  //! Class for computing approximate deformation gradient tensor
  class Compute_Deformation_Gradient : public PeridigmNS::Compute {

  public:
	
    //! Constructor.
    Compute_Deformation_Gradient( Teuchos::RCP<const Teuchos::ParameterList> params,
                             Teuchos::RCP<const Epetra_Comm> epetraComm_,
                             Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_);

    //! Destructor.
    ~Compute_Deformation_Gradient();

    //! Returns a vector of field IDs corresponding to the variables associated with the compute class.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Perform computation
    int compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks  ) const;

  private:

    // field ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_volumeFId;
    int m_horizonFId;
    int m_modelCoordinatesFId;
    int m_coordinatesFId;
    int m_shapeTensorInverseXXFId, m_shapeTensorInverseXYFId, m_shapeTensorInverseXZFId,
        m_shapeTensorInverseYXFId, m_shapeTensorInverseYYFId, m_shapeTensorInverseYZFId,
        m_shapeTensorInverseZXFId, m_shapeTensorInverseZYFId, m_shapeTensorInverseZZFId;
    int m_deformationGradientXXFId, m_deformationGradientXYFId, m_deformationGradientXZFId,
        m_deformationGradientYXFId, m_deformationGradientYYFId, m_deformationGradientYZFId,
        m_deformationGradientZXFId, m_deformationGradientZYFId, m_deformationGradientZZFId;

    // Unique ID for each instance of this class
    int myID;
    // Static member variable to generate IDs
    static int myIDGenerator;

  };

}

#endif // PERIDIGM_COMPUTE_DEFORMATION_GRADIENT_HPP
#endif // COMPUTE_CLASS
