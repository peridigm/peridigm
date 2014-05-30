/*! \file Peridigm_BoundaryCondition.cpp */

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

#include "Peridigm_BoundaryCondition.hpp"
#include "Peridigm.hpp"
#include <boost/math/special_functions/fpclassify.hpp>


using namespace std;

PeridigmNS::BoundaryCondition::BoundaryCondition(const string & name_,const Teuchos::ParameterList& bcParams_,Teuchos::RCP<Epetra_Vector> toVector_,Peridigm * peridigm_, const bool isCumulative_)
: peridigm(peridigm_),
  name(name_),
  toVector(toVector_),
  coord(0),
  muParserX(0.0),
  muParserY(0.0),
  muParserZ(0.0),
  muParserT(0.0),
  tensorOrder(SCALAR),
  isCumulative(isCumulative_)
{
  bcType = to_boundary_condition_type(bcParams_);
  string nodeSet = bcParams_.get<string>("Node Set");
  tidy_string(nodeSet);
  nodeSetName = nodeSet;
  coord = to_index(to_spatial_coordinate(bcParams_));
  function = bcParams_.get<string>("Value");

  // Set up muParser
  try {
    muParser.DefineVar("x", &muParserX);
    muParser.DefineVar("y", &muParserY);
    muParser.DefineVar("z", &muParserZ);
    muParser.DefineVar("t", &muParserT);
    muParser.DefineFun(_T("rnd"), mu::Rnd, false);
  }
  catch (mu::Parser::exception_type &e)
    TEUCHOS_TEST_FOR_EXCEPT_MSG(1, e.GetMsg());

  if(toVector->Map().ElementSize()==1)
  {
    tensorOrder = SCALAR;
    TEUCHOS_TEST_FOR_EXCEPTION(coord!=0,std::logic_error,"ERROR: The specified boundary condition coordinate must be X or omitted for BCs on scalar fields.");
  }
  else if(toVector->Map().ElementSize()==3)
    tensorOrder = VECTOR;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"ERROR: Boundary conditions have not been implemented for fields of tensor order.");
}

void PeridigmNS::BoundaryCondition::evaluateMuParser(const int & localNodeID, double & currentValue, double & previousValue, const double & timeCurrent, const double & timePrevious){
  Teuchos::RCP<Epetra_Vector> x = peridigm->getX();
  const Epetra_BlockMap& threeDimensionalMap = x->Map();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(threeDimensionalMap.ElementSize() != 3, "**** setVectorValues() must be called with map having element size = 3.\n");
  muParserX = (*x)[localNodeID*3];
  muParserY = (*x)[localNodeID*3 + 1];
  muParserZ = (*x)[localNodeID*3 + 2];
  muParserT = timeCurrent;
  try {
    currentValue = muParser.Eval();
  }
  catch (mu::Parser::exception_type &e)
  TEUCHOS_TEST_FOR_EXCEPT_MSG(1, e.GetMsg());

  muParserT = timePrevious;
  try {
    previousValue = muParser.Eval();
  }
  catch (mu::Parser::exception_type &e)
  TEUCHOS_TEST_FOR_EXCEPT_MSG(1, e.GetMsg());

  // if this is any other boundary condition besides prescribed displacement
  // get the previous value from evaluating the string function as above
  // if it is a perscribed displacement bc, check to see if the increment is
  // zero. This could happen if the user specifies a constant prescribed displacement.
  // If so, the increment should be the current prescribed value
  // minus the existing field value instead of the parser evaluation
  if(bcType==PRESCRIBED_DISPLACEMENT && currentValue - previousValue == 0.0)
  {
    Teuchos::RCP<Epetra_Vector> previousDisplacement = peridigm->getU();
    previousValue = (*previousDisplacement)[localNodeID*3 + coord];
  }

  // TODO: we should revisit how prescribed boundary conditions interact with initial conditions
  // in the case that the prescribed boundary condition at time zero doesn't match the inital condition
  // who wins?
}

PeridigmNS::DirichletBC::DirichletBC(const string & name_,const Teuchos::ParameterList& bcParams_,Teuchos::RCP<Epetra_Vector> toVector_,Peridigm * peridigm_,const bool isCumulative_)
: BoundaryCondition(name_,bcParams_,toVector_,peridigm_,isCumulative_){
}

void PeridigmNS::DirichletBC::apply(Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets, const double & timeCurrent, const double & timePrevious){
  // get the tensor order of the bc field:
  const int fieldDimension = to_dimension_size(tensorOrder);

  try{
    muParser.SetExpr(function);
  }
  catch (mu::Parser::exception_type &e)
  TEUCHOS_TEST_FOR_EXCEPT_MSG(1, e.GetMsg());

  // apply the bc to every element in the entire domain
  if(to_set_definition(nodeSetName)==FULL_DOMAIN)
  {
    for(int localNodeID = 0; localNodeID < toVector->MyLength(); localNodeID++) {
      double currentValue = 0.0;
      double previousValue = 0.0;
      evaluateMuParser(localNodeID,currentValue,previousValue,timeCurrent);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite(currentValue), "**** NaN returned by dirichlet BC evaluation.\n");
      if(isCumulative)
        (*toVector)[localNodeID*fieldDimension + coord] += currentValue;
      else
        (*toVector)[localNodeID*fieldDimension + coord] = currentValue;
    }
  }
  // apply the bc only to specific node sets
  else{
    std::map< std::string, std::vector<int> >::iterator itBegin;
    std::map< std::string, std::vector<int> >::iterator itEnd;
    if (to_set_definition(nodeSetName) == ALL_SETS){
      itBegin = nodeSets->begin();
      itEnd = nodeSets->end();
    }
    else{
      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(nodeSetName) == nodeSets->end(),
                                  "**** Error in DirichletBC::apply(), node set not found: " + nodeSetName + "\n");
      itBegin = nodeSets->find(nodeSetName);
      itEnd = itBegin; itEnd++;
    }
    for(std::map<std::string,std::vector<int> > ::iterator setIt=itBegin;setIt!=itEnd;++setIt){
      vector<int> & nodeList = setIt->second;
      for(unsigned int i=0 ; i<nodeList.size() ; i++){
        int localNodeID = toVector->Map().LID(nodeList[i]);
        if(localNodeID != -1) {
          double currentValue = 0.0;
          double previousValue = 0.0;
          evaluateMuParser(localNodeID,currentValue,previousValue,timeCurrent);
          TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite(currentValue), "**** NaN returned by dirichlet BC evaluation.\n");
          if(isCumulative)
            (*toVector)[localNodeID*fieldDimension + coord] += currentValue;
          else
            (*toVector)[localNodeID*fieldDimension + coord] = currentValue;
        }
      }
    }
  }
}

PeridigmNS::DirichletIncrementBC::DirichletIncrementBC(const string & name_,
  const Teuchos::ParameterList& bcParams_,
  Teuchos::RCP<Epetra_Vector> toVector_,
  Peridigm * peridigm_,
  const bool isCumulative_,
  const double & coeff_,
  const double & deltaTCoeff_)
: BoundaryCondition(name_,bcParams_,toVector_,peridigm_,isCumulative_),
  coeff(coeff_),
  deltaTCoeff(deltaTCoeff_){
}

void PeridigmNS::DirichletIncrementBC::apply(Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets, const double & timeCurrent, const double & timePrevious){
  // for temperature bcs, the previous time should always be zero
  const double timePrevious_ = bcType==PRESCRIBED_TEMPERATURE ? 0.0 : timePrevious;

  // get the tensor order of the bc field:
  const int fieldDimension = to_dimension_size(tensorOrder);

  try{
    muParser.SetExpr(function);
  }
  catch (mu::Parser::exception_type &e)
  TEUCHOS_TEST_FOR_EXCEPT_MSG(1, e.GetMsg());

  // apply the bc to every element in the entire domain
  if(to_set_definition(nodeSetName)==FULL_DOMAIN)
  {
    for(int localNodeID = 0; localNodeID < toVector->MyLength(); localNodeID++) {
      double currentValue = 0.0;
      double previousValue = 0.0;
      evaluateMuParser(localNodeID,currentValue,previousValue,timeCurrent,timePrevious_);
      const double value = coeff * (currentValue - previousValue)
                 + deltaTCoeff * (currentValue - previousValue) * (1.0 / (timeCurrent - timePrevious_));
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite(value), "**** NaN returned by dirichlet increment BC evaluation.\n");
      if(isCumulative)
        (*toVector)[localNodeID*fieldDimension + coord] += value;
      else
        (*toVector)[localNodeID*fieldDimension + coord] = value;
    }
  }
  else
  {
    std::map< std::string, std::vector<int> >::iterator itBegin;
    std::map< std::string, std::vector<int> >::iterator itEnd;
    if (to_set_definition(nodeSetName) == ALL_SETS){
      itBegin = nodeSets->begin();
      itEnd = nodeSets->end();
    }
    else{
      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(nodeSetName) == nodeSets->end(),
                                  "**** Error in DirichletIncrementBC::apply(), node set not found: " + nodeSetName + "\n");
      itBegin = nodeSets->find(nodeSetName);
      itEnd = itBegin; itEnd++;
    }
    for(std::map<std::string,std::vector<int> > ::iterator setIt=itBegin;setIt!=itEnd;++setIt){
      vector<int> & nodeList = setIt->second;
      for(unsigned int i=0 ; i<nodeList.size() ; i++){
        int localNodeID = toVector->Map().LID(nodeList[i]);
        if(localNodeID != -1) {
          double currentValue = 0.0;
          double previousValue = 0.0;
          evaluateMuParser(localNodeID,currentValue,previousValue,timeCurrent,timePrevious_);
          const double value = coeff * (currentValue - previousValue)
                         + deltaTCoeff * (currentValue - previousValue) * (1.0 / (timeCurrent - timePrevious_));
          TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite(value), "**** NaN returned by dirichlet increment BC evaluation.\n");
          if(isCumulative)
            (*toVector)[localNodeID*fieldDimension + coord] += value;
          else
            (*toVector)[localNodeID*fieldDimension + coord] = value;
        }
      }
    }
  }
}


