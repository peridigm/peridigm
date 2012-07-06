/*! \file Peridigm_ComputeManager.hpp */
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
#ifndef PERIDIGM_COMPUTEMANAGER_HPP
#define PERIDIGM_COMPUTEMANAGER_HPP

#include <vector>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Peridigm_Block.hpp>
#include <compute/Peridigm_Compute.hpp>

namespace PeridigmNS {
  
  class ComputeManager {
    
  public:
    
    //! Basic constructor.
    ComputeManager( Teuchos::RCP<Teuchos::ParameterList>& params, PeridigmNS::Peridigm *peridigm_ );
    
    //! Return list of field specs that the compute manager is handling
    std::vector<Field_NS::FieldSpec> getFieldSpecs();

    //! Destructor.
    virtual ~ComputeManager();

    //! Initialize the compute classes
    virtual void initialize(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks);

    //! Fire the individual compute objects
    virtual void compute(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks);

  private:
    
    //! Copy constructor.
    ComputeManager( const ComputeManager& CM );
    
    //! Assignment operator.
    ComputeManager& operator=( const ComputeManager& CM );

    //! Return valid input parameterlist
    Teuchos::ParameterList getValidParameterList();
    
    //! Individual compute objects
    std::vector< Teuchos::RCP<const PeridigmNS::Compute> > computeObjects;

    //! Parent pointer
    PeridigmNS::Peridigm *peridigm;

  };
  
}
 
#endif //PERIDIGM_COMPUTEMANAGER_HPP
