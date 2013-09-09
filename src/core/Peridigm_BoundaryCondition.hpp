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
#include "muParser/muParser.h"
#include "muParser/muParserPeridigmFunctions.h"
#include <Epetra_Vector.h>


using namespace std;

namespace PeridigmNS {

class Peridigm;

/*! \brief boundary and intial conditions (Initial conditions are boundary conditions).
 */
class BoundaryCondition {
public:

  //! Constructor.
  BoundaryCondition(const string & name_,const Teuchos::ParameterList& bcParams_,Teuchos::RCP<Epetra_Vector> bcVector_,Peridigm * peridigm_);

  //! Destructor.
  virtual ~BoundaryCondition(){}

  //! Give the name
  const string getName()const{return name;}

  //! give the name of the associated node set
  const string getNodeSetName()const{return nodeSetName;}

  //! give the type of boundary condition
  const Boundary_Condition_Type getType()const{return bcType;}

  //! give the field the bc is applied to
  Teuchos::RCP<Epetra_Vector> getBCVector()const{return bcVector;}

  //! give the coodinate of the bc
  const int getCoord()const{return coord;}

  //! apply the boundary condition
  virtual void apply(Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets, const double & timeCurrent=0.0, const double & timePrevious=0.0)=0;

  //! evaluate muParser
  void evaluateMuParser(const int & localNodeID, double & currentValue, double & previousValue, const double & timeCurrent=0.0, const double & timePrevious=0.0);

  //! Set the vector values
//  void setVectorValues_(Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets,
//    const double & timeCurrent=0.0,
//    const double & timePrevious=0.0,
//    const double & multiplier=1.0);

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
  Teuchos::RCP<Epetra_Vector> bcVector;

  //! coordinate direction associated with bc
  int coord;

  //! string defined funciton
  string function;

  //! Function parser
  mu::Parser muParser;

  //! @name Variables for function parser.
  //@{
  double muParserX;
  double muParserY;
  double muParserZ;
  double muParserT;
  //@}

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
  DirichletBC(const string & name_,const Teuchos::ParameterList& bcParams_,Teuchos::RCP<Epetra_Vector> bcVector_,Peridigm * peridigm_);

  //! Destructor.
  ~DirichletBC(){}

  //! apply the boundary condition
  virtual void apply(Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets, const double & timeCurrent=0.0, const double & timePrevious=0.0);

protected:

};

/*! \brief simply apply the function values to the vector
 */
class DirichletIncrementBC : public BoundaryCondition{
public:

  //! Constructor.
  DirichletIncrementBC(const string & name_,
    const Teuchos::ParameterList& bcParams_,
    Teuchos::RCP<Epetra_Vector> bcVector_,Peridigm * peridigm_,
    const double & coeff_=1.0,
    const double & deltaTCoeff_=0.0);

  //! Destructor.
  ~DirichletIncrementBC(){}

  //! apply the boundary condition
  virtual void apply(Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets, const double & timeCurrent=0.0, const double & timePrevious=0.0);

private:

  double coeff;
  double deltaTCoeff;
};






}

#endif // PERIDIGM_BOUNARYCONDITION_HPP
