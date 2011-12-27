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

// #include <Teuchos_RCP.hpp>
// #include <Epetra_Map.h>
// #include <Epetra_Vector.h>
// #include <Epetra_Import.h>
#include <Teuchos_ParameterList.hpp>

namespace PeridigmNS {

/*! \brief Processes boundary and intial conditions.
 *
 * Boundary conditions may be applied to any combination of node sets or blocks.
 *
 * Initial velocities
 * Initial displacements...is this really used?
 * Prescribed displacement boundary conditions (displacement as a function of time)
 * Body forces
 * Leave prescribed velocity boundary conditions (velocity as a function of time) for future work?
 *
 * Use functions for everything, at least from input deck point of view
 * Question:  How expensive is the expression parser?
 * Should create noop function as a poor-man's active periods
 * Allow for wacky IC/BC via compute-manager-style interface, a la user subroutines in Sierra/SM
 * Material state variables can already be initialized by a compute manager
 * Make sure compute manager initialize is called before material model initialize
 * Nix initial_conditions directory?
 * Move all muParser/randomNumber stuff out of Peridigm.cpp, into BC/IC classes, don't forget CMakeLists.txt
 * 
 * Initial velocity:  Applied to mothership vector once
 * Initial displacments:  Applied to mothership vector once
 * Prescribed displacement:  Applied to mothership vector, tangent, and displacement increment for QS/TD?
 * Body forces:  Applied to mothership vector at every time step
 *
 * Interface:
 * BoundaryAndInitialConditionManager(const Teuchos::RCP<Teuchos::ParameterList>& params)
 * initialize(Epetra_BlockMap)
 * applyInitialVelocities(Epetra_Vector vec)
 * applyInitialDisplacement(Epetra_Vector vec)
 * applyPrescribedDisplacement(Epetra_Vector vec)
 * applyPrescribedDisplacement(Epetra_CrsMatrix mat)
 * applyBodyForce(Epetra_Vector vec)
 */
  class BoundaryAndInitialConditionManager {
  public:

    //! Constructor
    BoundaryAndInitialConditionManager(const Teuchos::ParameterList& contactModelParams);

    //! Destructor
    ~BoundaryAndInitialConditionManager(){}

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
