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

#include <fstream>
#include <time.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <netcdf.h>
#include <exodusII.h>

#include <Epetra_Comm.h>
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include <Teuchos_Assert.hpp>

#include "Peridigm.hpp"
#include "Peridigm_SolverManager.hpp"

PeridigmNS::SolverManager::SolverManager(Teuchos::RCP<Teuchos::ParameterList> sParams, 
                                         PeridigmNS::Peridigm *peridigm_)
  : peridigm(peridigm_), sParams(sParams) {


  
  // No input to validate; no solver requested
  if (sParams == Teuchos::null) {
    return;
  }

  //TODO
  // Throws exception if parameters not present or of wrong type
  // Teuchos::ParameterList validator can't validate all input -- it mainly checks for presence of invalid input and invalid input types
  // Additional checking needed below 
  //Teuchos::ParameterList validParameterList = getValidParameterList();
  //bool isValid = true;
  //try {
    //params->validateParameters(validParameterList);
  //}
  //catch(Teuchos::Exceptions::InvalidParameterName &excpt)  {std::cout<<excpt.what(); isValid=false;}
  //catch(Teuchos::Exceptions::InvalidParameterType &excpt)  {std::cout<<excpt.what(); isValid=false;}
  //catch(Teuchos::Exceptions::InvalidParameterValue &excpt) {std::cout<<excpt.what(); isValid=false;}
  //catch(...) {isValid=false;}
  //if (!isValid){
    //std::cout.flush();
    //TEUCHOS_TEST_FOR_EXCEPTION(1, std::invalid_argument, "PeridigmNS::SolverManager:::SolverManager() -- Invalid parameter, type or value.");
  //}

}

Teuchos::ParameterList PeridigmNS::SolverManager::getValidParameterList() {

  // prevent Teuchos from converting parameter types
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes intParam(false), dblParam(false), strParam(false);
  intParam.allowInt(true);
  dblParam.allowDouble(true);
  strParam.allowString(true);

  // Construct a ParameterList containing valid entries for Solver
  Teuchos::ParameterList validParameterList("Solver");
  setDoubleParameter("Initial Time",0.0,"Start time for this Solver",&validParameterList,dblParam);
  setDoubleParameter("Final Time",0.0,"Finish time for this Solver",&validParameterList,dblParam);
  validParameterList.set("Verbose",false);

  Teuchos::ParameterList& validSolverVerletParameterList = validParameterList.sublist("Verlet");
  setDoubleParameter("Safety Factor",1.0,"Time step safety factor",&validSolverVerletParameterList,dblParam);
  setDoubleParameter("Fixed dt",1.0,"User defined fixed time step",&validSolverVerletParameterList,dblParam);

  Teuchos::ParameterList& validSolverQSParameterList = validParameterList.sublist("QuasiStatic");

  return validParameterList;
}

PeridigmNS::SolverManager::~SolverManager() {
}

void PeridigmNS::SolverManager::executeSolver() {

    peridigm->execute(sParams);
}

