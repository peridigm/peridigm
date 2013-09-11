/*! \file Peridigm_BoundaryAndInitialConditionManager.hpp */

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

#ifndef PERIDIGM_BOUNARYANDINITIALCONDITIONMANAGER_HPP
#define PERIDIGM_BOUNARYANDINITIALCONDITIONMANAGER_HPP

#include <Epetra_Vector.h>
#include <Epetra_FECrsMatrix.h>
#include <Teuchos_ParameterList.hpp>

#include "Peridigm_Discretization.hpp"
#include "Peridigm_BoundaryCondition.hpp"

#include <vector>

using namespace std;

namespace PeridigmNS {

  class Peridigm;

/*! \brief Processes boundary and intial conditions.
 */
  class BoundaryAndInitialConditionManager {
  public:

    //! Constructor.
    BoundaryAndInitialConditionManager(const Teuchos::ParameterList& boundaryAndInitialConditionParams,Peridigm * parent);

    //! Destructor.
    ~BoundaryAndInitialConditionManager(){}

    //! Initialize boundary conditions, etc.
    void initialize(Teuchos::RCP<Discretization> discretization);

    //! Initialize the node sets on the bc manager
    void initializeNodeSets(Teuchos::RCP<Discretization> discretization);

    //! Get node sets.
    Teuchos::RCP< std::map< std::string, std::vector<int> > > getNodeSets() {
      return nodeSets;
    }

    //! Apply all initial conditions
    void applyInitialConditions();

    //! Apply all boundary conditions
    void applyBoundaryConditions(const double & timeCurrent=0.0, const double & timePrevious=0.0);

    //! Apply all boundary conditions
    void applyForceContributions(const double & timeCurrent=0.0, const double & timePrevious=0.0);

    //! Apply all boundary conditions
    void clearForceContributions();

    //! Update the current coordinates
    void updateCurrentCoordinates();

    //! Copies entries corresponding to kinematic boundary contitions into the vector of reaction forces.
    void applyKinematicBC_ComputeReactions(Teuchos::RCP<const Epetra_Vector> force, Teuchos::RCP<Epetra_Vector> reaction);

    //! Set rows corresponding to kinematic boundary conditions to zero.
    void applyKinematicBC_InsertZeros(Teuchos::RCP<Epetra_Vector> vec);

    //! Set rows and columns corresponding to kinematic boundary conditions to zero and put 1.0 on the diagonal.
    void applyKinematicBC_InsertZerosAndSetDiagonal(Teuchos::RCP<Epetra_FECrsMatrix> mat);

  protected:

    //! Boundary and initial condition parameters
    Teuchos::ParameterList params;

    //! Node sets
    Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets;

    //! Determine if a string is the name of a file; if so return the file name, if not return an empty string.
    std::string nodeSetStringToFileName(std::string str);

    //! Ref your parent instantiator
    Peridigm  * peridigm;

    //! Set of all the boundary conditions
    vector<Teuchos::RCP<BoundaryCondition> > boundaryConditions;

    //! Set of all the initial conditions
    vector<Teuchos::RCP<BoundaryCondition> > initialConditions;

    //! Set of all the force contributions
    vector<Teuchos::RCP<BoundaryCondition> > forceContributions;



  private:

    // Private to prohibit use.
    BoundaryAndInitialConditionManager();

    // Private to prohibit use.
    BoundaryAndInitialConditionManager(const BoundaryAndInitialConditionManager&);

    // Private to prohibit use.
    BoundaryAndInitialConditionManager& operator=(const BoundaryAndInitialConditionManager&);

  };

}

#endif // PERIDIGM_BOUNARYANDINITIALCONDITIONMANAGER_HPP
