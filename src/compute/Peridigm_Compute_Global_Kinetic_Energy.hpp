/*! \file Peridigm_Compute_Global_Kinetic_Energy.hpp */

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

ComputeClass(Global_Kinetic_Energy,Compute_Global_Kinetic_Energy)

#else

#ifndef PERIDIGM_COMPUTE_GLOBAL_KINETIC_ENERGY_HPP
#define PERIDIGM_COMPUTE_GLOBAL_KINETIC_ENERGY_HPP

#include "Peridigm_Compute.hpp"
#include "Peridigm_Compute_Kinetic_Energy.hpp"

namespace PeridigmNS {

  //! Class for filling the kinetic energy vector
  class Compute_Global_Kinetic_Energy : public PeridigmNS::Compute_Kinetic_Energy {

  public:

    //! Constructor.
    Compute_Global_Kinetic_Energy( Teuchos::RCP<const Teuchos::ParameterList> params,
                             Teuchos::RCP<const Epetra_Comm> epetraComm_,
                             Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_);

    //! Destructor.
    ~Compute_Global_Kinetic_Energy();

    //! Perform computation
    virtual int compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks  ) const;

    //! Returns a vector of field IDs corresponding to the variables associated with the compute class.
    virtual std::vector<int> FieldIds() const { return Compute_Kinetic_Energy::FieldIds(); }
  };
}

#endif // PERIDIGM_COMPUTE_GLOBAL_KINETIC_ENERGY_HPP
#endif // COMPUTE_CLASS
