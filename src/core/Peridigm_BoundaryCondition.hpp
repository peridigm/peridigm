/*! \file Peridigm_BoundaryCondition.hpp */

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

#ifndef PERIDIGM_BOUNARYCONDITION_HPP
#define PERIDIGM_BOUNARYCONDITION_HPP

#include "Peridigm_Enums.hpp"
#include <Epetra_Vector.h>

#include <Trilinos_version.h>
#if TRILINOS_MAJOR_MINOR_VERSION >= 111100
#include "RTC_FunctionRTC.hh"
#else
#include "FunctionRTC.hh"
#endif

using namespace std;

namespace PeridigmNS {

class Peridigm;

/*! \brief boundary and intial conditions (Initial conditions are boundary conditions).
 */
class BoundaryCondition {
public:

  //! Constructor.
  BoundaryCondition(const string & name_,
                    const Teuchos::ParameterList& bcParams_,
                    Teuchos::RCP<Epetra_Vector> toVector_,
                    Peridigm * peridigm_);

  //! Destructor.
  virtual ~BoundaryCondition(){}

  //! Give the name
  string getName()const{return name;}

  //! give the name of the associated node set
  string getNodeSetName()const{return nodeSetName;}

  //! give the type of boundary condition
  Boundary_Condition_Type getType()const{return bcType;}

  //! give the field the bc is applied to
  Teuchos::RCP<Epetra_Vector> getBCVector()const{return toVector;}

  //! give the coodinate of the bc
  int getCoord()const{return coord;}

  //! apply the boundary condition
  virtual void apply(Teuchos::RCP< std::map< std::string,
                     std::vector<int> > > nodeSets,
                     const double & timeCurrent = 0.0,
                     const double & timePrevious = 0.0) = 0;

  //! evaluate function parser
  void evaluateParser(const int & localNodeID,
                      double & currentValue,
                      double & previousValue,
                      const double & timeCurrent = 0.0,
                      const double & timePrevious = 0.0);

protected:

  //! Ref your parent instantiator
  Peridigm  * peridigm;

  //! string name
  string name;

  //! Node set name
  string nodeSetName;

  //! BC Type
  Boundary_Condition_Type bcType;

  //! the field to apply the boundary condtion to
  Teuchos::RCP<Epetra_Vector> toVector;

  //! coordinate direction associated with bc
  int coord;

  //! string defined funciton
  string function;

  //! Run-time compiler, used as function parser
  Teuchos::RCP<PG_RuntimeCompiler::Function> rtcFunction;

  Tensor_Order tensorOrder;

private:

  // Private to prohibit use.
  BoundaryCondition();

  // Private to prohibit use.
  BoundaryCondition(const BoundaryCondition&);

  // Private to prohibit use.
  BoundaryCondition& operator=(const BoundaryCondition&);
};


/*! \brief simply apply the function values to the vector
 */
class DirichletBC : public BoundaryCondition{
public:

  //! Constructor.
  DirichletBC(const string & name_,
              const Teuchos::ParameterList& bcParams_,
              Teuchos::RCP<Epetra_Vector> toVector_,
              Peridigm * peridigm_);

  //! Destructor.
  ~DirichletBC(){}

  //! apply the boundary condition
  virtual void apply(Teuchos::RCP< std::map< std::string,
                     std::vector<int> > > nodeSets,
                     const double & timeCurrent = 0.0,
                     const double & timePrevious = 0.0);
};

/*! \brief simply apply the function values to the vector
 */
class DirichletIncrementBC : public BoundaryCondition{
public:

  //! Constructor.
  DirichletIncrementBC(const string & name_,
                       const Teuchos::ParameterList& bcParams_,
                       Teuchos::RCP<Epetra_Vector> toVector_,
                       Peridigm * peridigm_,
                       const double & coeff_ = 1.0,
                       const double & deltaTCoeff_ = 0.0,
                       bool computeChangeRelativeToInitialValue_ = false);

  //! Destructor.
  ~DirichletIncrementBC(){}

  //! apply the boundary condition
  virtual void apply(Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets,
                     const double & timeCurrent = 0.0,
                     const double & timePrevious = 0.0);

private:

  double coeff;
  double deltaTCoeff;
  bool computeChangeRelativeToInitialValue;
};

}

#endif // PERIDIGM_BOUNARYCONDITION_HPP
