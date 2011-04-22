//*! \file Peridigm_OutputManager.cpp */
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

#include <time.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <Epetra_Comm.h>
#include <Teuchos_TestForException.hpp>
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Peridigm_OutputManager_VTK_XML.hpp"
#include <Field.h>
#include "PdVTK.h"

PeridigmNS::OutputManager_VTK_XML::OutputManager_VTK_XML(const Teuchos::RCP<Teuchos::ParameterList>& params) {

  // No input to validate; no output requested
  if (params == Teuchos::null) {
    iWrite = false;
    return;
  }

  // Throws exception if parameters not present or of wrong type
  // Teuchos::ParameterList validator can't validate all input -- it mainly checks for presence of invalid input and invalid input types
  // Additional checking needed below 
  Teuchos::ParameterList validParameterList = getValidParameterList();
  bool isValid = true;
  try {
    params->validateParameters(validParameterList);
  }
  catch(Teuchos::Exceptions::InvalidParameterName &excpt)  {std::cout<<excpt.what(); isValid=false;}
  catch(Teuchos::Exceptions::InvalidParameterType &excpt)  {std::cout<<excpt.what(); isValid=false;}
  catch(Teuchos::Exceptions::InvalidParameterValue &excpt) {std::cout<<excpt.what(); isValid=false;}
  catch(...) {isValid=false;}
  if (!isValid) TEST_FOR_EXCEPTION(1, std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- Invalid parameter, type or value.");

  try {
    numProc = params->INVALID_TEMPLATE_QUALIFIER get<int>("NumProc");
  }
  catch ( const std::exception& e) {
    TEST_FOR_EXCEPTION(1, std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- numProc not present.");
  }

  try {
    myPID = params->INVALID_TEMPLATE_QUALIFIER get<int>("MyPID");
  }
  catch ( const std::exception& e) {
    TEST_FOR_EXCEPTION(1,  std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- MyPID not present.");
  }

  // Default to no output
  frequency = params->get<int>("Output Frequency",-1); 

  // Default to BINARY output
  outputFormat = params->get<string>("Output Format","BINARY"); 
  TEST_FOR_EXCEPTION( (outputFormat != "ASCII") && (outputFormat != "BINARY"),  std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- Unknown output format. Must be ASCII or BINARY.");

  // Default to not write full neighborlist
  writeNeighborlist = params->get<bool>("Bond Family",false); 
  TEST_FOR_EXCEPTION( (numProc != 1) && (writeNeighborlist),  std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- Parallel write of bond families not currently supported.");

  // Output filename base
  filenameBase = params->get<string>("Output Filename","dump"); 

  // User-requested fields for output 
  materialOutputFields = sublist(params, "Material Output Fields");

  // Initialize count (number of times write() has been called)
  // Initialize to -1 because first call to write() corresponds to timstep 0
  count = -1;

  // With VTK, every object writes
  iWrite = true;

  // VTK doesn't like spaces or . so replace them with underscore
  {
  int warningFlag = 0;
  string outString;
  outString.append("\n\n***WARNING***\n");
  outString.append("PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- Avoid use of filenames containing '.' (period) and ' ' (space) with VTK.\n"); 
  outString.append("Changing ");
  outString.append(filenameBase);
  outString.append(" to ");
  for ( unsigned int i = 0; i < filenameBase.length(); i++) {
    if (filenameBase[i] ==' ' || filenameBase[i]=='.')  {
      filenameBase.replace(i,1,"_");
      warningFlag = 1;
    }
  }
  outString.append(filenameBase);
  outString.append(".\n\n\n");
  if (warningFlag) std::cout << outString; 
  }

  // Create VTK collection writer
  if (outputFormat == "ASCII")
    vtkWriter = Teuchos::rcp(new PdVTK::CollectionWriter(filenameBase.c_str(), numProc, myPID, PdVTK::vtkASCII));
  else if (outputFormat == "BINARY")
    vtkWriter = Teuchos::rcp(new PdVTK::CollectionWriter(filenameBase.c_str(), numProc, myPID, PdVTK::vtkBINARY));

}

Teuchos::ParameterList PeridigmNS::OutputManager_VTK_XML::getValidParameterList() {

  // prevent Teuchos from converting parameter types
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes intParam(false), dblParam(false), strParam(false);
  intParam.allowInt(true);
  dblParam.allowDouble(true);
  strParam.allowString(true);

  Teuchos::ParameterList validParameterList("Output");
  setIntParameter("MyPID",0,"Process ID",&validParameterList,intParam);
  setIntParameter("NumProc",0,"Number of Process IDs",&validParameterList,intParam);
  validParameterList.set("Output File Type","VTK_XML");
  validParameterList.set("Output Filename","dump");
  Teuchos::setStringToIntegralParameter<int>("Output Format","BINARY","ASCII or BINARY",Teuchos::tuple<string>("ASCII","BINARY"),&validParameterList);
  setIntParameter("Output Frequency",-1,"Frequency of Output",&validParameterList,intParam);
  validParameterList.set("Parallel Write",true);
  Teuchos::ParameterList& matFields = validParameterList.sublist("Material Output Fields");
  { // Valid output fields for Linear Elastic material type
  Teuchos::ParameterList& matType = matFields.sublist("Linear Elastic");
  matType.set("Displacement",true);
  matType.set("Velocity",true);
  matType.set("Acceleration",true);
  matType.set("Force Density",true);
  matType.set("Contact Force Density",true);
  matType.set("Dilatation",true);
  matType.set("ID",true);
  matType.set("Proc Num",true);
  matType.set("Weighted Volume",true);
  matType.set("Damage",true);
  matType.set("Volume",true);
  }
  { // Valid output fields for Elastic Plastic material type
  Teuchos::ParameterList& matType = matFields.sublist("Elastic Plastic");
  matType.set("Displacement",true);
  matType.set("Velocity",true);
  matType.set("Acceleration",true);
  matType.set("Force Density",true);
  matType.set("Contact Force Density",true);
  matType.set("Dilatation",true);
  matType.set("ID",true);
  matType.set("Proc Num",true);
  matType.set("Weighted Volume",true);
  matType.set("Damage",true);
  matType.set("Lambda",true);
  matType.set("Volume",true);
  matType.set("Shear_Correction_Factor",true);
  }

  return validParameterList;
}

PeridigmNS::OutputManager_VTK_XML::~OutputManager_VTK_XML() {

  if (!iWrite) return;

  // get current system time
  time_t rawtime;
  struct tm *timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  // Construct comment string to be written to master .pvd file
  string comment;
  comment.append("Peridigm Version XXX: "); 
  comment.append(asctime(timeinfo));

  vtkWriter->close(comment);  
}

void PeridigmNS::OutputManager_VTK_XML::write(Teuchos::RCP<const Epetra_Vector> x,
                                              Teuchos::RCP<const Epetra_Vector> u,
                                              Teuchos::RCP<const Epetra_Vector> v,
                                              Teuchos::RCP<const Epetra_Vector> a,
                                              Teuchos::RCP<const Epetra_Vector> force,
                                              Teuchos::RCP<PeridigmNS::DataManager> dataManager,
                                              Teuchos::RCP<const NeighborhoodData> neighborhoodData,
                                              Teuchos::RCP<Teuchos::ParameterList>& forceStateDesc) {
  if (!iWrite) return;

  // increment index count
  count = count + 1;

  // Only write if frequency count match
  if (frequency<=0 || count%frequency!=0) return;

  // Initialize grid if needed
  static int rebalanceCount = 0;
  if (grid.GetPointer() == NULL || rebalanceCount != dataManager->getRebalanceCount()) {
    double *xptr;
    x->ExtractView( &xptr );
    int length = (x->Map()).NumMyElements();
    grid = PdVTK::getGrid(xptr,length);
    rebalanceCount = dataManager->getRebalanceCount();
  }

  // Currenly only support "Linear Elastic" and "Elastic Plastic" material models
  Teuchos::RCP<Teuchos::ParameterList> thisMaterial;
  if (materialOutputFields->isSublist ("Linear Elastic")) {
    thisMaterial = sublist(materialOutputFields, "Linear Elastic");
  }
  else if (materialOutputFields->isSublist ("Elastic Plastic")) {
    thisMaterial = sublist(materialOutputFields, "Elastic Plastic");
  }
  else // Unrecognized material model. Throw error.
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
    "Peridigm::OutputManager: Unknown material model. Only \"Linear Elastic\" or \"Elastic Plastic\" currently supported.");

  if (thisMaterial->isParameter("Volume")) {
	  double *volume;
	  dataManager->getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&volume);
	  PdVTK::writeField<double>(grid,Field_NS::VOLUME,volume);
  }

  if (thisMaterial->isParameter("Displacement")) {
	  double *uptr;
	  u->ExtractView( &uptr );
	  PdVTK::writeField<double>(grid,Field_NS::DISPL3D,uptr);
  }

  if (thisMaterial->isParameter("Velocity")) {
   double *vptr;
   v->ExtractView( &vptr );
   PdVTK::writeField<double>(grid,Field_NS::VELOC3D,vptr);
  }

  if (thisMaterial->isParameter("Acceleration")) {
   double *aptr;
   a->ExtractView( &aptr );
   PdVTK::writeField<double>(grid,Field_NS::ACCEL3D,aptr);
  }

  if (thisMaterial->isParameter("Force Density")) {
    double *fptr;
    force->ExtractView( &fptr );
    PdVTK::writeField<double>(grid,Field_NS::FORCE_DENSITY3D,fptr);
  }

  if (thisMaterial->isParameter("ID")) {
    // Get map corresponding to x
    const Epetra_BlockMap& xMap = x->Map();
    PdVTK::writeField<int>(grid,Field_NS::ID,xMap.MyGlobalElements());
  }

  std::vector<int> proc_num;
  if (thisMaterial->isParameter("Proc Num")) {
    // Get map corresponding to x
    const Epetra_BlockMap& xMap = x->Map();
    int length = xMap.NumMyElements();
    proc_num.assign (length,myPID); 
    PdVTK::writeField<int>(grid,Field_NS::PROC_NUM,&proc_num[0]);
  }

  // Currenly only support "Linear Elastic" and "Elastic Plastic" material models
  Teuchos::RCP<Teuchos::ParameterList> thisForceState;
  if (materialOutputFields->isSublist ("Linear Elastic")) {
    thisForceState = sublist(forceStateDesc, "Linear Elastic");
  } 
  else if (materialOutputFields->isSublist ("Elastic Plastic")) {
    thisForceState = sublist(forceStateDesc, "Elastic Plastic");
  }
  else // Unrecognized material model. Throw error.
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
    "Peridigm::OutputManager: Unknown material model. Only \"Linear Elastic\" or \"Elastic Plastic\" currently supported.");

  Teuchos::ParameterList::ConstIterator iter;
  string outstring;
  for(iter = thisForceState->begin(); iter != thisForceState->end(); iter++) {
    if (thisForceState->isParameter(thisMaterial->name(iter))) {
      outstring = thisMaterial->name(iter);

      double *dataptr;
      if (outstring == "Dilatation") {
        dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dataptr);
        PdVTK::writeField<double>(grid,Field_NS::DILATATION,dataptr);
      }
      else if (outstring == "Weighted_Volume") {
        dataManager->getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&dataptr);
        PdVTK::writeField<double>(grid,Field_NS::WEIGHTED_VOLUME,dataptr);
      }
      else if (outstring == "Damage") {
        dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dataptr);
        PdVTK::writeField<double>(grid,Field_NS::DAMAGE,dataptr);
      }
      else if (outstring == "Lambda") {
        dataManager->getData(Field_NS::LAMBDA, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dataptr);
        PdVTK::writeField<double>(grid,Field_NS::LAMBDA,dataptr);
      }
      else {
        // Unknown field
      }

    }
  }

  // All pointers reset; now write data
  double current_time = forceStateDesc->get<double>("Time");
  vtkWriter->writeTimeStep(current_time,grid);
//  vtkWriter->writeTimeStep(count,grid);
}

void PeridigmNS::OutputManager_VTK_XML::write(Teuchos::RCP<PeridigmNS::DataManager> dataManager,
                                              Teuchos::RCP<const NeighborhoodData> neighborhoodData,
                                              Teuchos::RCP<Teuchos::ParameterList>& forceStateDesc) {
  if (!iWrite) return;

  // increment index count
  count = count + 1;

  // Only write if frequency count match
  if (frequency<=0 || count%frequency!=0) return;

  // Initialize grid if needed
  static int rebalanceCount = 0;
  if (grid.GetPointer() == NULL || rebalanceCount != dataManager->getRebalanceCount()) {
    double *xptr;
    Teuchos::RCP<Epetra_Vector> myX =  dataManager->getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE);
    myX->ExtractView( &xptr );
    //int length = (myX->Map()).NumMyElements();
    // Use only the number of owned elements
    int length = (dataManager->getOwnedIDVectorMap())->NumMyElements();
    grid = PdVTK::getGrid(xptr,length);
    rebalanceCount = dataManager->getRebalanceCount();
  }

  // Currenly only support "Linear Elastic" and "Elastic Plastic" material models
  Teuchos::RCP<Teuchos::ParameterList> thisMaterial;
  if (materialOutputFields->isSublist ("Linear Elastic")) {
    thisMaterial = sublist(materialOutputFields, "Linear Elastic");
  }
  else if (materialOutputFields->isSublist ("Elastic Plastic")) {
    thisMaterial = sublist(materialOutputFields, "Elastic Plastic");
  }
  else // Unrecognized material model. Throw error.
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
    "Peridigm::OutputManager: Unknown material model. Only \"Linear Elastic\" or \"Elastic Plastic\" currently supported.");

  if (thisMaterial->isParameter("Volume")) {
    double *volume;
    dataManager->getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&volume);
    PdVTK::writeField<double>(grid,Field_NS::VOLUME,volume);
  }

  if (thisMaterial->isParameter("Displacement")) {
    double *uptr;
    dataManager->getData(Field_NS::DISPL3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView( &uptr );
    PdVTK::writeField<double>(grid,Field_NS::DISPL3D,uptr);
  }

  if (thisMaterial->isParameter("Velocity")) {
    double *vptr;
    dataManager->getData(Field_NS::VELOC3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView( &vptr );
    PdVTK::writeField<double>(grid,Field_NS::VELOC3D,vptr);
  }

  if (thisMaterial->isParameter("Acceleration")) {
// MLP: This is currently not stored in datamanager
    double *aptr;
    dataManager->getData(Field_NS::ACCEL3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView( &aptr );
    PdVTK::writeField<double>(grid,Field_NS::ACCEL3D,aptr);
  }

  if (thisMaterial->isParameter("Force Density")) {
    double *fptr;
    dataManager->getData(Field_NS::FORCE_DENSITY3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView( &fptr );
    PdVTK::writeField<double>(grid,Field_NS::FORCE_DENSITY3D,fptr);
  }

  if (thisMaterial->isParameter("ID")) {
    // Get map corresponding to x
    Teuchos::RCP<Epetra_Vector> myX =  dataManager->getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE);
    const Epetra_BlockMap& xMap = myX->Map();
    PdVTK::writeField<int>(grid,Field_NS::ID,xMap.MyGlobalElements());
  }

  std::vector<int> proc_num;
  if (thisMaterial->isParameter("Proc Num")) {
    // Get map corresponding to x
    Teuchos::RCP<Epetra_Vector> myX =  dataManager->getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE);
//    const Epetra_BlockMap& xMap = myX->Map();
//    int length = xMap.NumMyElements();
    // Use only the number of owned elements
    int length = (dataManager->getOwnedIDVectorMap())->NumMyElements();
    proc_num.assign (length,myPID);
    PdVTK::writeField<int>(grid,Field_NS::PROC_NUM,&proc_num[0]);
  }

  // Currenly only support "Linear Elastic" and "Elastic Plastic" material models
  Teuchos::RCP<Teuchos::ParameterList> thisForceState;
  if (materialOutputFields->isSublist ("Linear Elastic")) {
    thisForceState = sublist(forceStateDesc, "Linear Elastic");
  }
  else if (materialOutputFields->isSublist ("Elastic Plastic")) {
    thisForceState = sublist(forceStateDesc, "Elastic Plastic");
  }
  else // Unrecognized material model. Throw error.
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
    "Peridigm::OutputManager: Unknown material model. Only \"Linear Elastic\" or \"Elastic Plastic\" currently supported.");

  Teuchos::ParameterList::ConstIterator iter;
  string outstring;
  for(iter = thisForceState->begin(); iter != thisForceState->end(); iter++) {
    if (thisForceState->isParameter(thisMaterial->name(iter))) {
      outstring = thisMaterial->name(iter);

      double *dataptr;
      if (outstring == "Dilatation") {
        dataManager->getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dataptr);
        PdVTK::writeField<double>(grid,Field_NS::DILATATION,dataptr);
      }
      else if (outstring == "Weighted_Volume") {
        dataManager->getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&dataptr);
        PdVTK::writeField<double>(grid,Field_NS::WEIGHTED_VOLUME,dataptr);
      }
      else if (outstring == "Damage") {
        dataManager->getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dataptr);
        PdVTK::writeField<double>(grid,Field_NS::DAMAGE,dataptr);
      }
      else if (outstring == "Lambda") {
        dataManager->getData(Field_NS::LAMBDA, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dataptr);
        PdVTK::writeField<double>(grid,Field_NS::LAMBDA,dataptr);
      }
      else {
        // Unknown field
      }

    }
  }

  // All pointers reset; now write data
  double current_time = forceStateDesc->get<double>("Time");
  vtkWriter->writeTimeStep(current_time,grid);
//  vtkWriter->writeTimeStep(count,grid);
}
