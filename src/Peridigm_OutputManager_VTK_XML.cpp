//*! \file Peridigm_OutputManager.cpp */
// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

#include <time.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <Epetra_Comm.h>
#include <Teuchos_TestForException.hpp>

#include "Peridigm_OutputManager_VTK_XML.hpp"
#include <Field.h>

PeridigmNS::OutputManager_VTK_XML::OutputManager_VTK_XML(const Teuchos::RCP<Teuchos::ParameterList>& params) {
  // MLP: Insert call to parameterlist validator here

  // Throws exception if parameters not present
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
    TEST_FOR_EXCEPTION(1,  std::invalid_argument, "PeridigmNS::OutputManager_VTK_XML:::OutputManager_VTK_XML() -- myPID not present.");
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

PeridigmNS::OutputManager_VTK_XML::~OutputManager_VTK_XML() {

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
                                              Teuchos::RCP<const Epetra_MultiVector> scalarConstitutiveData,
                                              Teuchos::RCP<const NeighborhoodData> neighborhoodData,
                                              Teuchos::RCP<Teuchos::ParameterList>& forceStateDesc) {
  if (!iWrite) return;

  // increment index count
  count = count + 1;

  // Only write if frequency count match
  if (frequency<=0 || count%frequency!=0) return;

  // Initialize grid if needed
  if (grid.GetPointer() == NULL) {
    double *xptr;
    x->ExtractView( &xptr );
    int length = (x->Map()).NumMyElements();
    grid = PdVTK::getGrid(xptr,length);
  }

  // MLP: For now, hardcode to "Linear Elastic" material type. Changes to facilitate additional material types will be refelected here.
  Teuchos::RCP<Teuchos::ParameterList> thisMaterial = sublist(materialOutputFields, "Linear Elastic");

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

  // MLP: For now, hardcode to "Linear Elastic" material type. Changes in Peridigm_ModelEvaluator to facilitate additional material types will be 
  // refelected here.
  Teuchos::RCP<Teuchos::ParameterList> thisForceState = sublist(forceStateDesc, "Linear Elastic");

  Teuchos::ParameterList::ConstIterator iter;
  string outstring;
  for(iter = thisMaterial->begin(); iter != thisMaterial->end(); iter++) {
    if (thisForceState->isParameter(thisMaterial->name(iter))) {
      outstring = thisMaterial->name(iter);
      // get column index into scalarConstitutiveData for this parameter
      int colIdx = Teuchos::getValue<int>(thisForceState->getEntry(thisMaterial->name(iter)));

      double *dataptr;
      dataptr = (*scalarConstitutiveData)[colIdx];
      if (outstring == "Dilatation") {
        PdVTK::writeField<double>(grid,Field_NS::DILATATION,dataptr);
      }
      else if (outstring == "Weighted Volume") {
        PdVTK::writeField<double>(grid,Field_NS::WEIGHTED_VOLUME,dataptr);
      }
      else if (outstring == "Damage") {
        PdVTK::writeField<double>(grid,Field_NS::DAMAGE,dataptr);
      }
      else {
        // Unknown field
      }

    }
  }

  // All pointers reset; now write data
  vtkWriter->writeTimeStep(count,grid);

}
