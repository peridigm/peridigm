/*! \file Peridigm_OutputManager_Container.hpp */
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
#ifndef PERIDIGM_OUTPUTMANAGER_CONTAINER_HPP
#define PERIDIGM_OUTPUTMANAGER_CONTAINER_HPP

#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Peridigm_Block.hpp>
#include <Peridigm_OutputManager.hpp>

namespace PeridigmNS {

  class OutputManagerContainer {

  public:

    //! Basic constructor.
    OutputManagerContainer(){};

    //! Destructor.
    ~OutputManagerContainer() {};

    //! Add new output manager
    void add(Teuchos::RCP<PeridigmNS::OutputManager> OM) {
      outputManagers.push_back( OM );
    }

    //! Write to all output managers in container
    void write(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks, double current_time) {
      std::vector< Teuchos::RCP< PeridigmNS::OutputManager > >::iterator it;
      for ( it=outputManagers.begin() ; it < outputManagers.end(); it++ )
        (*it)->write(blocks, current_time);
    }

    //! Multiply output frequency of all output managers in container
    //  for the sake of reducing load step size in Adaptive Quasi-static
    void multiplyOutputFrequency(double multiplier){
      std::vector< Teuchos::RCP< PeridigmNS::OutputManager > >::iterator it;
      for ( it=outputManagers.begin() ; it < outputManagers.end(); it++ )
        (*it)->multiplyOutputFrequency(multiplier);
    }

    //! Change output frequency of all output managers in container
    //  for the sake of switch from Quasi-static to explicit solver
    void changeOutputFrequency(int output_frequency){
      std::vector< Teuchos::RCP< PeridigmNS::OutputManager > >::iterator it;
      for ( it=outputManagers.begin() ; it < outputManagers.end(); it++ )
        (*it)->changeOutputFrequency(output_frequency);
    }

  protected:

    //! Container for RCPs to individual output managers
    std::vector< Teuchos::RCP< PeridigmNS::OutputManager > > outputManagers;

  private:

    //! Copy constructor.
    OutputManagerContainer( const OutputManagerContainer& OMC );

    //! Assignment operator.
    OutputManagerContainer& operator=( const OutputManagerContainer& OMC );

  };

}

#endif //PERIDIGM_OUTPUTMANAGER_CONTAINER_HPP
