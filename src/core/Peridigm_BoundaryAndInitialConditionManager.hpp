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

#include "muParser/muParserDef.h" // \todo Is this required?
#include "muParser/muParser.h"

#include "Peridigm_AbstractDiscretization.hpp"

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
 * Should we have separate managers for boundary and initial conditions?
 * Change "all node sets" to "all blocks"
 * How do we really want to handle nonzero initial displacements? 
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

    //! Constructor.
    BoundaryAndInitialConditionManager(const Teuchos::ParameterList& boundaryAndInitialConditionParams);

    //! Destructor.
    ~BoundaryAndInitialConditionManager(){}

    //! Initialize node sets, etc.
    void initialize(Teuchos::RCP<AbstractDiscretization> discretization);

    //! Apply initial displacements.
    void applyInitialDisplacements(Teuchos::RCP<Epetra_Vector> x,
                                   Teuchos::RCP<Epetra_Vector> u,
                                   Teuchos::RCP<Epetra_Vector> y);

    //! Apply initial velocities.
    void applyInitialVelocities(Teuchos::RCP<const Epetra_Vector> x,
                                Teuchos::RCP<Epetra_Vector> v);

    /** \brief Set velocity in rows corresponding to kinematic boundary conditions.
     ** Intended for use with velocity vector in explicit time integration such that
     ** the explicit integration scheme will propagate these values through to the
     ** displacement and current coordinate vectors. **/
    void applyKinematicBC_SetVelocity(double timeCurrent,
                                      double timePrevious,
                                      Teuchos::RCP<const Epetra_Vector> x,
                                      Teuchos::RCP<Epetra_Vector> vec);

    /** \brief Set displacement in rows corresponding to kinematic boundary conditions.
     ** Intended use is setting nonzero displacements at time zero. **/
    void applyKinematicBC_SetDisplacement(double timeCurrent,
                                          Teuchos::RCP<const Epetra_Vector> x,
                                          Teuchos::RCP<Epetra_Vector> vec);


    /** \brief Set displacement increment in rows corresponding to kinematic boundary conditions.
     ** Intended for use with displacement increment vector in quasi-static time integration. **/
    void applyKinematicBC_SetDisplacementIncrement(double timeCurrent,
                                                   double timePrevious,
                                                   Teuchos::RCP<const Epetra_Vector> x,
                                                   Teuchos::RCP<Epetra_Vector> vec);

    //! Set rows corresponding to kinematic boundary conditions to zero.
    void applyKinematicBC_InsertZeros(Teuchos::RCP<Epetra_Vector> vec);

    //! Set rows and columns corresponding to kinematic boundary conditions to zero and put 1.0 on the diagonal.
    void applyKinematicBC_InsertZerosAndPutOnesOnDiagonal(Teuchos::RCP<Epetra_FECrsMatrix> mat);

  protected:

    //! Boundary and initial condition parameters
    Teuchos::ParameterList params;

    // \todo Define nodesets with an Epetra_BlockMap, then export to x, v, etc., and rebalance as needed

    //! Node sets
    Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets;

    //! Function parser
    mu::Parser muParser;

    //! @name Variables for function parser.
    //@{ 
    double muParserX;
    double muParserY;
    double muParserZ;
    double muParserT;
    //@}

    //! Set either a displacement or a displacement increment.
    void setVectorValues(double timeCurrent,
                         double timePrevious,
                         Teuchos::RCP<const Epetra_Vector> x,
                         Teuchos::RCP<Epetra_Vector> vec,
                         bool setIncrement,
                         double multiplier);

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
