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
#include <cmath>
#include <sstream>


using namespace std;

PeridigmNS::BoundaryCondition::BoundaryCondition(const string & name_,
                                                 const Teuchos::ParameterList& bcParams_,
                                                 Teuchos::RCP<Epetra_Vector> toVector_,
                                                 Peridigm * peridigm_)
: peridigm(peridigm_),
  name(name_),
  toVector(toVector_),
  coord(0),
  tensorOrder(SCALAR)
{
  Teuchos::ParameterList bcParams(bcParams_);
  bcType = to_boundary_condition_type(bcParams);
  string nodeSet = bcParams.get<string>("Node Set");
  tidy_string(nodeSet);
  nodeSetName = nodeSet;
  coord = to_index(to_spatial_coordinate(bcParams));
  if(bcParams.isType<double>("Value")){
    std::stringstream ss;
    ss << bcParams.get<double>("Value");
    bcParams.set("Value", ss.str());
  }
  function = bcParams.get<string>("Value");

  // set up RTCompiler
  rtcFunction = Teuchos::rcp<PG_RuntimeCompiler::Function>(new PG_RuntimeCompiler::Function(5, "rtcBoundaryConditionFunction"));
  rtcFunction->addVar("double", "x");
  rtcFunction->addVar("double", "y");
  rtcFunction->addVar("double", "z");
  rtcFunction->addVar("double", "t");
  rtcFunction->addVar("double", "value");

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

void PeridigmNS::BoundaryCondition::evaluateParser(const int & localNodeID, double & currentValue, double & previousValue, const double & timeCurrent, const double & timePrevious){
  Teuchos::RCP<Epetra_Vector> x = peridigm->getX();
  const Epetra_BlockMap& threeDimensionalMap = x->Map();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(threeDimensionalMap.ElementSize() != 3, "**** setVectorValues() must be called with map having element size = 3.\n");

  bool success(true);
  // set the coordinates and set the return value to 0.0
  if(success)
    rtcFunction->varValueFill(0, (*x)[localNodeID*3]);
  if(success)
    success = rtcFunction->varValueFill(1, (*x)[localNodeID*3 + 1]);
  if(success)
    success = rtcFunction->varValueFill(2, (*x)[localNodeID*3 + 2]);
  if(success)
    success = rtcFunction->varValueFill(4, 0.0);
  // evaluate at previous time
  if(success)
    success = rtcFunction->varValueFill(3, timePrevious);
  if(success)
    success = rtcFunction->execute();
  if(success)
    previousValue = rtcFunction->getValueOfVar("value");
  // evaluate at current time
  if(success)
    success = rtcFunction->varValueFill(3, timeCurrent);
  if(success)
    success = rtcFunction->execute();
  if(success)
    currentValue = rtcFunction->getValueOfVar("value");
  if(!success){
    string msg = "\n**** Error in BoundaryCondition::evaluateParser().\n";
    msg += "**** " + rtcFunction->getErrors() + "\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!success, msg);
  }

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

PeridigmNS::DirichletBC::DirichletBC(const string & name_,
                                     const Teuchos::ParameterList& bcParams_,
                                     Teuchos::RCP<Epetra_Vector> toVector_,
                                     Peridigm * peridigm_)
: BoundaryCondition(name_,bcParams_,toVector_,peridigm_){
}

void PeridigmNS::DirichletBC::apply(Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets,
                                    const double & timeCurrent,
                                    const double & timePrevious){

  // get the tensor order of the bc field:
  const int fieldDimension = to_dimension_size(tensorOrder);

  string rtcFunctionString = function;
  if(rtcFunctionString.find("value") == string::npos)
    rtcFunctionString = "value = " + rtcFunctionString;
  bool success = rtcFunction->addBody(rtcFunctionString);
  if(!success){
    string msg = "\n**** Error:  rtcFunction->addBody(function) returned error code in PeridigmNS::DirichletBC::apply().\n";
    msg += "**** " + rtcFunction->getErrors() + "\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!success, msg);
  }

  // apply the bc to every element in the entire domain
  if(to_set_definition(nodeSetName)==FULL_DOMAIN)
  {
    for(int localNodeID = 0; localNodeID < toVector->Map().NumMyElements(); localNodeID++) {
      double currentValue = 0.0;
      double previousValue = 0.0;
      evaluateParser(localNodeID,currentValue,previousValue,timeCurrent);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(currentValue), "**** NaN returned by dirichlet BC evaluation.\n");
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
          evaluateParser(localNodeID,currentValue,previousValue,timeCurrent);
          TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(currentValue), "**** NaN returned by dirichlet BC evaluation.\n");
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
  const double & coeff_,
  const double & deltaTCoeff_,
  bool computeChangeRelativeToInitialValue_)
: BoundaryCondition(name_,bcParams_,toVector_,peridigm_),
  coeff(coeff_),
  deltaTCoeff(deltaTCoeff_),
  computeChangeRelativeToInitialValue(computeChangeRelativeToInitialValue_){
}

void PeridigmNS::DirichletIncrementBC::apply(Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets,
                                             const double & timeCurrent,
                                             const double & timePrevious){

  // if setting BC at the onset of a simulation, timeCurrent will equal timePrevious
  // do not apply incremental bc in this case
  if(timeCurrent == timePrevious)
    return;

  const double timePrevious_ = computeChangeRelativeToInitialValue ? 0.0 : timePrevious;

  // get the tensor order of the bc field:
  const int fieldDimension = to_dimension_size(tensorOrder);

  string rtcFunctionString = function;
  if(rtcFunctionString.find("value") == string::npos)
    rtcFunctionString = "value = " + rtcFunctionString;
  bool success = rtcFunction->addBody(rtcFunctionString);
  if(!success){
    string msg = "\n**** Error:  rtcFunction->addBody(function) returned nonzero error code in PeridigmNS::DirichletBC::apply().\n";
    msg += "**** " + rtcFunction->getErrors() + "\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!success, msg);
  }

  // apply the bc to every element in the entire domain
  if(to_set_definition(nodeSetName)==FULL_DOMAIN)
  {
    for(int localNodeID = 0; localNodeID < toVector->Map().NumMyElements(); localNodeID++) {
      double currentValue = 0.0;
      double previousValue = 0.0;
      evaluateParser(localNodeID,currentValue,previousValue,timeCurrent,timePrevious_);
      const double value = coeff * (currentValue - previousValue)
                 + deltaTCoeff * (currentValue - previousValue) * (1.0 / (timeCurrent - timePrevious_));
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(value), "**** NaN returned by dirichlet increment BC evaluation.\n");
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
          evaluateParser(localNodeID,currentValue,previousValue,timeCurrent,timePrevious_);
          const double value = coeff * (currentValue - previousValue)
                         + deltaTCoeff * (currentValue - previousValue) * (1.0 / (timeCurrent - timePrevious_));
          TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(value), "**** NaN returned by dirichlet increment BC evaluation.\n");
          (*toVector)[localNodeID*fieldDimension + coord] = value;
        }
      }
    }
  }
}


