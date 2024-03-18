/*! \file Peridigm_HorizonManager.cpp */

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

#include "Peridigm_HorizonManager.hpp"
#include <Teuchos_Assert.hpp>
#include <Teuchos_Exceptions.hpp>
#include <iterator>

using std::vector;
using std::map;
using std::pair;
using std::copy;
using std::cout;
using std::string;
using std::ofstream;
using std::istringstream;
using std::stringstream;
using std::endl;
using std::setprecision;
using std::istream_iterator;
using std::back_inserter;


PeridigmNS::HorizonManager::HorizonManager() {
  // set up RTCompiler
  rtcFunction = Teuchos::rcp<PG_RuntimeCompiler::Function>(new PG_RuntimeCompiler::Function(4, "rtcHorizonFunction"));
  rtcFunction->addVar("double", "x");
  rtcFunction->addVar("double", "y");
  rtcFunction->addVar("double", "z");
  rtcFunction->addVar("double", "value");
}

PeridigmNS::HorizonManager& PeridigmNS::HorizonManager::self() {
  static HorizonManager horizonManager;
  return horizonManager;
}

void PeridigmNS::HorizonManager::loadHorizonInformationFromBlockParameters(Teuchos::ParameterList& blockParams) {

  // Find the horizon value for each block and record the default horizon value (if any)
  for(Teuchos::ParameterList::ConstIterator it = blockParams.begin() ; it != blockParams.end() ; it++){
    Teuchos::ParameterList& params = blockParams.sublist(it->first);
    bool hasConstantHorizon = params.isType<double>("Horizon");
    bool hasVariableHorizon = params.isType<string>("Horizon");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(hasConstantHorizon && hasVariableHorizon, "\n**** Error parsing horizon information!  Multiple horizon definitions found!\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!hasConstantHorizon && !hasVariableHorizon, "\n**** Error parsing horizon information!  Either the horizon is not defined or defined incorrectly!\n");

    // Record the horizon as a string regardless of whether or not it is constant
    string horizonString;
    if(hasConstantHorizon){
      double constantHorizon = params.get<double>("Horizon");
      stringstream horizonStringStream;
      horizonStringStream.precision(16);
      horizonStringStream << constantHorizon;
      horizonString = horizonStringStream.str();
    }
    else{
      horizonString = params.get<string>("Horizon");
    }

    // Parse space-delimited list of block names
    string blockNamesString = params.get<string>("Block Names");
    istringstream iss(blockNamesString);
    vector<string> blockNames;
    copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter<vector<string> >(blockNames));

    for(vector<string>::const_iterator it = blockNames.begin() ; it != blockNames.end() ; ++it){
      if( *it == "Default" || *it == "default" || *it == "DEFAULT" ){
        horizonIsConstant["default"] = hasConstantHorizon;
        horizonStrings["default"] = horizonString;
      }
      else{
        horizonIsConstant[*it] = hasConstantHorizon;
        horizonStrings[*it] = horizonString;
      }
    }
  }
}

bool PeridigmNS::HorizonManager::blockHasConstantHorizon(string blockName){
  bool isConstant;
  if(horizonIsConstant.find(blockName) != horizonIsConstant.end())
    isConstant = horizonIsConstant[blockName];
  else if(horizonIsConstant.find("default") != horizonIsConstant.end())
    isConstant = horizonIsConstant["default"];
  else{
    string msg = "\n**** Error, no Horizon parameter found for block " + blockName + " and no default block parameter list provided.\n";
    TEUCHOS_TEST_FOR_TERMINATION(true, msg);
  }
  return isConstant;
}

double PeridigmNS::HorizonManager::getBlockConstantHorizonValue(string blockName){
  string name;
  if(horizonStrings.find(blockName) != horizonStrings.end())
    name = blockName;
  else if(horizonStrings.find("default") != horizonStrings.end())
    name = "default";
  else{
    string msg = "\n**** Error, no Horizon parameter found for block " + blockName + " and no default block parameter list provided.\n";
    TEUCHOS_TEST_FOR_TERMINATION(true, msg);
  }
  double x(0.0), y(0.0), z(0.0);
  double horizon = evaluateHorizon(name, x, y, z);
  return horizon;
}

double PeridigmNS::HorizonManager::evaluateHorizon(string blockName, double x, double y, double z){
  string name;
  if(horizonStrings.find(blockName) != horizonStrings.end())
    name = blockName;
  else if(horizonStrings.find("default") != horizonStrings.end())
    name = "default";
  else{
    string msg = "\n**** Error, no Horizon parameter found for block " + blockName + " and no default block parameter list provided.\n";
    TEUCHOS_TEST_FOR_TERMINATION(true, msg);
  }
  string horizonFunction = horizonStrings[name];
  double horizonValue(0.0);

  string rtcFunctionString = horizonFunction;
  if(rtcFunctionString.find("value") == string::npos)
    rtcFunctionString = "value = " + rtcFunctionString;
  bool success = rtcFunction->addBody(rtcFunctionString);
  if(success)
    rtcFunction->varValueFill(0, x);
  if(success)
    success = rtcFunction->varValueFill(1, y);
  if(success)
    success = rtcFunction->varValueFill(2, z);
  if(success)
    success = rtcFunction->varValueFill(3, 0.0);
  if(success)
    success = rtcFunction->execute();
  if(success)
    horizonValue = rtcFunction->getValueOfVar("value");
  if(!success){
    string msg = "\n**** Error in HorizonManager::evaluateHorizon().\n";
    msg += "**** " + rtcFunction->getErrors() + "\n";
    TEUCHOS_TEST_FOR_TERMINATION(!success, msg);
  }

  return horizonValue;
}




// void PeridigmNS::Peridigm::checkHorizon(Teuchos::RCP<Discretization> peridigmDisc, map<string, double> & horizonStrings){
//   // Warn the user if it appears the horizon is too large and may cause memory to fill
//   const double warningPerc = 0.80; // TODO: This might need to be adjusted
//   const int mesh_size = peridigmDisc->getNumElem();
//   const bool mesh_size_large = mesh_size > 10000; // TODO This might need to be adjusted
//   const double maxPerc = (double)peridigmDisc->getMaxNumBondsPerElem() / (double)mesh_size;
//   if(maxPerc >= warningPerc && mesh_size_large){
//     if(peridigmComm->MyPID() == 0){
//       cout << "** Warning: elements were detected with large neighborhood sizes relative to the number of local elements.\n"
//            << "** The largest element neighborhood contains " << maxPerc * 100 << "% of the elements on this processor.\n"
//            << "** This may indicate the horizon was selected too large and could lead to memory capacity being exceeded.\n\n";
//     }
//     return;
//   }
//   // Warn the user if a block's horizon is too large compared to the max element diameter
//   const double warningSizeRatio = 25.0; // TODO: This might need to be adjusted
//   string blockName = "";
//   double horizonValue = 0.0;
//   const double maxRad = peridigmDisc->getMaxElementRadius();
//   for(map<string, double>::const_iterator it = horizonStrings.begin() ; it != horizonStrings.end() ; it++){
//     if(it->second > warningSizeRatio * maxRad && mesh_size_large){
//       blockName = it->first;
//       horizonValue = it->second;
//       if(peridigmComm->MyPID() == 0){
//         cout << "** Warning: The horizon for " << blockName << " is " << horizonValue << ", which is\n"
//              << "** more than " << warningSizeRatio <<  " times the max element radius (" << maxRad << ").\n"
//              << "** This may indicate the horizon was selected too large\n"
//              << "** and could lead to memory capacity being exceeded.\n\n";
//       }
//     }
//   }
// }
