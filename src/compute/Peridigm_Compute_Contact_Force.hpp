/*! \file Peridigm_Compute_Contact_Force.hpp */

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

ComputeClass(Contact_Force,Compute_Contact_Force)

#else

#ifndef PERIDIGM_COMPUTE_CONTACT_FORCE_HPP
#define PERIDIGM_COMPUTE_CONTACT_FORCE_HPP

#include "Peridigm_Compute.hpp"

namespace PeridigmNS {

  //! Class for filling acceleration vector
  class Compute_Contact_Force : public PeridigmNS::Compute {

  public:
	
    //! Constructor.
    Compute_Contact_Force( Teuchos::RCP<const Epetra_Comm> epetraComm_ );

    //! Destructor.
    ~Compute_Contact_Force();

    //! Returns a vector of field IDs corresponding to the variables associated with the compute class.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Perform computation
    int compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks  ) const;

  private:

    // field ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_volumeFieldId;
    int m_contactForceDensityFieldId;
    int m_contactForceFieldId;
  };
}

#endif // PERIDIGM_COMPUTE_CONTACT_FORCE_HPP
#endif // COMPUTE_CLASS
