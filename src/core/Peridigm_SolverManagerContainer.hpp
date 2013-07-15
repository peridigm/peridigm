/*! \file Peridigm_SolverManager_Container.hpp */
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
#ifndef PERIDIGM_SOLVERMANAGER_CONTAINER_HPP
#define PERIDIGM_SOLVERMANAGER_CONTAINER_HPP

#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Peridigm_SolverManager.hpp>

namespace PeridigmNS {

  class SolverManagerContainer {

  public:

    //! Basic constructor.
    SolverManagerContainer(){};

    //! Destructor.
    ~SolverManagerContainer() {};

    //! Add new output manager
    void add(Teuchos::RCP<PeridigmNS::SolverManager> SM) {
      solverManagers.push_back( SM );
    }


    //! Evaluate solver for all solver managers in container
    void executeSolvers(){
      std::vector< Teuchos::RCP<PeridigmNS::SolverManager> >::iterator it;
      for ( it=solverManagers.begin() ; it < solverManagers.end(); it++ ){
        (*it)->executeSolver();
      }
    }

  protected:

    //! Container for RCPs to individual output managers
    std::vector< Teuchos::RCP<PeridigmNS::SolverManager> > solverManagers;

  private:

    //! Copy constructor.
    SolverManagerContainer( const SolverManagerContainer& SMC );

    //! Assignment operator.
    SolverManagerContainer& operator=( const SolverManagerContainer& SMC );

  };

}

#endif //PERIDIGM_SOLVERMANAGER_CONTAINER_HPP
