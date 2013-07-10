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
#include "muParser/muParser.h"

namespace PeridigmNS {

/*! \brief Processes boundary and intial conditions.
 */
  class BoundaryAndInitialConditionManager {
  public:

    //! Constructor.
    BoundaryAndInitialConditionManager(const Teuchos::ParameterList& boundaryAndInitialConditionParams);

    //! Destructor.
    ~BoundaryAndInitialConditionManager(){}

    //! Initialize node sets, etc.
    void initialize(Teuchos::RCP<Discretization> discretization);

    //! Get node sets.
    Teuchos::RCP< std::map< std::string, std::vector<int> > > getNodeSets() {
      return nodeSets;
    }

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

    /** \brief Determines temperature change at every node based on a user-defined temperature field. **/
    void applyTemperatureChange(double timeCurrent,
                                Teuchos::RCP<const Epetra_Vector> x,
                                Teuchos::RCP<Epetra_Vector> deltaT);

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

    //! Function parser
    mu::Parser muParser;

    //! @name Variables for function parser.
    //@{ 
    double muParserX;
    double muParserY;
    double muParserZ;
    double muParserT;
    //@}

    //! Flag indicating presence of thermal field.
    bool m_hasThermal;

    //! Set either a displacement or a displacement increment.
    void setVectorValues(double timeCurrent,
                         double timePrevious,
                         Teuchos::RCP<const Epetra_Vector> x,
                         Teuchos::RCP<Epetra_Vector> vec,
                         bool setIncrement,
                         double multiplier);

    //! Determine if a string is the name of a file; if so return the file name, if not return an empty string.
    std::string nodeSetStringToFileName(std::string str);

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
