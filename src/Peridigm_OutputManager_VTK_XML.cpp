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

#include <fstream>
#include <time.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <Epetra_Comm.h>
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include <Teuchos_TestForException.hpp>

#include "Peridigm_OutputManager_VTK_XML.hpp"
#include "mesh_output/vtk/Field.h"
#include "mesh_output/vtk/PdVTK.h"
#include "Peridigm.hpp"

PeridigmNS::OutputManager_VTK_XML::OutputManager_VTK_XML(const Teuchos::RCP<Teuchos::ParameterList>& params, 
                                                         PeridigmNS::Peridigm *peridigm_,
                                                         Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks) {
 
  // Assign parent pointer
  peridigm = peridigm_;
 
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

  // Create container to hold processor ID # data
  proc_num = Teuchos::rcp( new std::vector<int>() );

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

  // Flag if more than one block
  if (blocks->size() > 1)
    isMultiBlock = true;
  else 
    isMultiBlock = false;

  // Create vector of VTK collection writers, one per block
  vtkWriters.resize(blocks->size());
  // Iterators over vector of grids and over blocks. Each vector is the same size by above line.
  std::vector<PeridigmNS::Block>::iterator blockIt;
  std::vector<Teuchos::RCP<PdVTK::CollectionWriter> >::iterator writerIt;
  if (outputFormat == "ASCII") {
    for(blockIt = blocks->begin(), writerIt=vtkWriters.begin() ; blockIt != blocks->end() ; blockIt++, writerIt++)
      (*writerIt) = Teuchos::rcp(new PdVTK::CollectionWriter(filenameBase.c_str(), numProc, myPID, PdVTK::vtkASCII));
  }
  else if (outputFormat == "BINARY") {
    for(blockIt = blocks->begin(), writerIt=vtkWriters.begin() ; blockIt != blocks->end() ; blockIt++, writerIt++)
      (*writerIt) = Teuchos::rcp(new PdVTK::CollectionWriter(filenameBase.c_str(), numProc, myPID, PdVTK::vtkBINARY));
  }

}

Teuchos::ParameterList PeridigmNS::OutputManager_VTK_XML::getValidParameterList() {

  //! Todo: This code assumes knowledage of materials in material library. Replace this code when material model manager in place.

  // prevent Teuchos from converting parameter types
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes intParam(false), dblParam(false), strParam(false);
  intParam.allowInt(true);
  dblParam.allowDouble(true);
  strParam.allowString(true);

  // Get valid output fields from parent (Peridigm object)
  std::vector<Field_NS::FieldSpec> peridigmSpecs = peridigm->getFieldSpecs();

  // Get valid output fields from compute manager
  std::vector<Field_NS::FieldSpec> computeSpecs = peridigm->computeManager->getFieldSpecs();

  // Container for valid output fields from material classes
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > materialSpecs;

  // Construct a valid output parameter list based upon instantiated marterial and compute objects
  Teuchos::ParameterList validParameterList("Output");
  setIntParameter("MyPID",0,"Process ID",&validParameterList,intParam);
  setIntParameter("NumProc",0,"Number of Process IDs",&validParameterList,intParam);
  validParameterList.set("Output File Type","VTK_XML");
  validParameterList.set("Output Filename","dump");
  Teuchos::setStringToIntegralParameter<int>("Output Format","BINARY","ASCII or BINARY",Teuchos::tuple<string>("ASCII","BINARY"),&validParameterList);
  setIntParameter("Output Frequency",-1,"Frequency of Output",&validParameterList,intParam);
  validParameterList.set("Parallel Write",true);
  Teuchos::ParameterList& matFields = validParameterList.sublist("Material Output Fields");

  // Now loop over all instantiated materials, filling material output fields sublist for each material type
  for(unsigned int i=0; i<peridigm->materialModels->size() ; ++i){
    // Get name of underlying material
    string name = peridigm->materialModels->operator[](i)->Name();
    // Create sublist with that name
    Teuchos::ParameterList& matType = matFields.sublist(name);
    // Container to hold all the valid specs for this material type
    std::vector<Field_NS::FieldSpec> matTypeSpecs;
    // Aggregate specs from Peridigm object, ComputeManager object, and Material object
    // Do not insert any fieldSpecs with FieldLength == BOND (e.g., no bond data)
    for(unsigned int j=0; j<peridigmSpecs.size() ; ++j) {
      if (peridigmSpecs[j].getRelation() != Field_ENUM::BOND)
        matTypeSpecs.insert(matTypeSpecs.end(),peridigmSpecs[j]);
    }
    for(unsigned int j=0; j<computeSpecs.size() ; ++j) {
      if (computeSpecs[j].getRelation() != Field_ENUM::BOND)
        matTypeSpecs.insert(matTypeSpecs.end(),computeSpecs[j]);
    }
    materialSpecs = peridigm->materialModels->operator[](i)->VariableSpecs();
    for(unsigned int j=0; j<materialSpecs->size() ; ++j) {
      if (materialSpecs->operator[](j).getRelation() != Field_ENUM::BOND)
        matTypeSpecs.insert(matTypeSpecs.end(),materialSpecs->operator[](j));
    }
    // ID and ProcNum can be determined from any *Petra vector, so list them as well
    matTypeSpecs.insert(matTypeSpecs.end(),Field_NS::GID);
    matTypeSpecs.insert(matTypeSpecs.end(),Field_NS::PROC_NUM);
    // Remove duplicates
    std::unique(matTypeSpecs.begin(), matTypeSpecs.end());
    // Sort for consistency
    std::sort(matTypeSpecs.begin(), matTypeSpecs.end());
    // Now walk the matTypeSpec vector and create the parameterlist
    for(unsigned int j=0; j<matTypeSpecs.size() ; ++j) matType.set(matTypeSpecs.operator[](j).getLabel(),true);
    // Clear vector
    matTypeSpecs.clear();
  }

  return validParameterList;
}

PeridigmNS::OutputManager_VTK_XML::~OutputManager_VTK_XML() {

  if (!iWrite) return;

  // Only proc 0 writes the .pvd file
  if (myPID != 0) return;

  // get current system time
  time_t rawtime;
  struct tm *timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  // Construct comment string to be written to master .pvd file
  string comment;
  comment.append("Peridigm Version XXX: "); 
  comment.append(asctime(timeinfo));

  if (!isMultiBlock) // write normal .pvd file
    vtkWriters.front()->close(comment);
  else { // write pvd file as a container of vtm files
    string outFileName(filenameBase);
    outFileName += ".pvd";

    std::fstream outFile;
    outFile.open(outFileName.c_str(), fstream::out);

    outFile << "<?xml version=\"1.0\"?>" << std::endl;
    outFile << "<!--" << std::endl;
    outFile << comment << std::endl;
    outFile << "-->" << std::endl;
    outFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    outFile << "  <Collection>" << std::endl;

    deque<double> times = vtkWriters.front()->getTimes();
    for(unsigned int i=0;i<times.size();i++) {
      std::stringstream ss;
      ss << filenameBase.c_str() << "_t" << i << ".vtm";
      string vtmFileName(ss.str());
      outFile << "    " << "<DataSet timestep=\"" << times[i] << "\" " << "file=\"" << vtmFileName << "\"/>" << std::endl;
    }
    outFile << "  </Collection>" << std::endl;
    outFile << "</VTKFile>" << std::endl;

    outFile.close();

  } // end of else 

}

void PeridigmNS::OutputManager_VTK_XML::write(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks, double current_time) {

  if (!iWrite) return;

  // increment index count
  count = count + 1;

  // Only write if frequency count match
  if (frequency<=0 || count%frequency!=0) return;

  // Initialize/reinitialize grids if needed
  // Each block is always rebalanced at the same time, so each datamanager should always return the same
  // rebalance count. Hence, we keep only a single static int for the rebalance count.
  static int rebalanceCount = 0;
  bool reInitGrids = false;
  // If the first block rebalanced since last write, then all of them did. Force reinit of all grids.
  if (rebalanceCount != blocks->begin()->getDataManager()->getRebalanceCount()) {
    reInitGrids = true;
    rebalanceCount = blocks->begin()->getDataManager()->getRebalanceCount();
  }
  // Be sure vector of grids sized correctly
  if (blocks->size() != grids.size() ) {
    grids.resize(blocks->size());
    reInitGrids = true;
  }
  // Iterators over vector of grids and over blocks. Each vector is  same size by above line.
  std::vector<PeridigmNS::Block>::iterator blockIt;
  std::vector<vtkSmartPointer<vtkUnstructuredGrid> >::iterator gridIt;
  for(blockIt = blocks->begin(), gridIt=grids.begin() ; blockIt != blocks->end() ; blockIt++, gridIt++) {
    if (gridIt->GetPointer() == NULL || reInitGrids) {
      double *xptr;
      Teuchos::RCP<Epetra_Vector> myX =  blockIt->getDataManager()->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE);
      myX->ExtractView( &xptr );
      // Use only the number of owned elements
      int length = (blockIt->getDataManager()->getOwnedScalarPointMap())->NumMyElements();
      (*gridIt) = PdVTK::getGrid(xptr,length);
    }
  }

  if (blocks->size() == 1) // if only one block, don't write with block IDs
    this->write(blocks->begin()->getDataManager(),blocks->begin()->getNeighborhoodData(),grids.front(),vtkWriters.front(),current_time);
  else { // there are many blocks; write output using block IDs
    std::vector<Teuchos::RCP<PdVTK::CollectionWriter> >::iterator writerIt;
    for(blockIt = blocks->begin(), gridIt=grids.begin(), writerIt=vtkWriters.begin() ; blockIt != blocks->end() ; blockIt++, gridIt++, writerIt++)
      this->write(blockIt->getDataManager(),blockIt->getNeighborhoodData(),*gridIt,*writerIt,current_time,blockIt->getID());
  }

  // If multiblock simulation, write container file for each block
  if ( (blocks->size() > 1) && (myPID==0) ) {
    writeVTMFile(blocks);
  }

}

void PeridigmNS::OutputManager_VTK_XML::write(Teuchos::RCP<PeridigmNS::DataManager> dataManager,
                                              Teuchos::RCP<const NeighborhoodData> neighborhoodData,
                                              vtkSmartPointer<vtkUnstructuredGrid> grid,
                                              Teuchos::RCP<PdVTK::CollectionWriter> vtkWriter,
                                              double current_time, int block_id) {

  // Call compute manager; Updated any computed quantities before write
  // \todo This will break for multiple materials.
  peridigm->computeManager->compute(dataManager);

  Teuchos::ParameterList::ConstIterator i1;
  // Loop over the material types in the materialOutputFields parameterlist
  for (i1 = materialOutputFields->begin(); i1 != materialOutputFields->end(); ++i1) {
    const Teuchos::ParameterEntry& val1 = materialOutputFields->entry(i1);
    // const std::string& name1 = materialOutputFields->name(i1);
    // For each material type, loop over requested output fields and hook up pointers
    if (val1.isList()) { // each material type is a sublist
      const Teuchos::ParameterList& sublist = Teuchos::getValue<Teuchos::ParameterList>(val1);
      Teuchos::ParameterList::ConstIterator i2;
      for (i2=sublist.begin(); i2 != sublist.end(); ++i2) {
        const std::string& nm = sublist.name(i2);
        // use field name to get reference to const fieldSpec
        std::map<string, Field_NS::FieldSpec>::const_iterator i3; 
        i3 = Field_NS::FieldSpecMap::Map.find(nm); // Can't use operator[] on a const std::map
        TEST_FOR_EXCEPT_MSG(i3 == Field_NS::FieldSpecMap::Map.end(), "Failed to find reference to fieldSpec!");
        Field_NS::FieldSpec const &fs = i3->second;
        double *ptr; ptr = NULL;
        if (fs == Field_NS::GID) { // Handle special case of ID (int type)
          // Get map corresponding to x (COORD3D FieldSpec guaranteed to exist by Peridigm object)
          Teuchos::RCP<Epetra_Vector> myX = dataManager->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE);
          const Epetra_BlockMap& xMap = myX->Map();
          // hook up pointer to data
          PdVTK::writeField<int>(grid,Field_NS::GID,xMap.MyGlobalElements());
        }
         else if (fs == Field_NS::PROC_NUM) { // Handle special case of Proc_Num (int type)
          // Get map corresponding to x (COORD3D FieldSpec guaranteed to exist by Peridigm object)
          Teuchos::RCP<Epetra_Vector> myX =  dataManager->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE);
          // Use only the number of owned elements
          int length = (dataManager->getOwnedScalarPointMap())->NumMyElements();
          // If the length is zero, this means there are no on-processor points for this block
          if(length > 0)
            proc_num->assign(length,myPID);
          else
            proc_num->assign(1,myPID); // Avoids access error in subsequent call to proc_num->at(0)
          // hook up pointer to data
          PdVTK::writeField<int>(grid,Field_NS::PROC_NUM,&(proc_num->at(0)));
       }
        else { // Handle all other cases (double type)
          if (fs.get_temporal() != Field_ENUM::TWO_STEP) // If stateless, get STEP_NONE
            dataManager->getData(fs, Field_ENUM::STEP_NONE)->ExtractView(&ptr);
          else // If stateful, get STEP_NP1
            dataManager->getData(fs, Field_ENUM::STEP_NP1)->ExtractView(&ptr);
          // hook up pointer to data
          PdVTK::writeField<double>(grid,fs,ptr);
        }
      }
    }
  }

  // All pointers reset; now write data
  vtkWriter->writeTimeStep(current_time,grid,block_id);
}


void PeridigmNS::OutputManager_VTK_XML::writeVTMFile(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks) {
  std::stringstream ss;
  std::size_t index = vtkWriters.front()->getIndex();
  ss << filenameBase.c_str() << "_t" << index << ".vtm";
  string outFileName(ss.str());

  std::fstream outFile;
  outFile.open(outFileName.c_str(), fstream::out); 

  outFile << "<?xml version=\"1.0\"?>" << std::endl;
  outFile << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
  outFile << "  <vtkMultiBlockDataSet>" << std::endl;

  std::vector<PeridigmNS::Block>::iterator blockIt;
  std::vector<Teuchos::RCP<PdVTK::CollectionWriter> >::iterator writerIt;
  for(blockIt = blocks->begin(), writerIt=vtkWriters.begin() ; blockIt != blocks->end() ; blockIt++, writerIt++) {
    outFile << "     <Block index=\"" << blockIt->getID() << "\">" << std::endl;
    for(int i=0;i<numProc;i++) {
      std::stringstream ss;
      ss << filenameBase << "_b" << blockIt->getID() << "_t" << (*writerIt)->getIndex() << "_" << i << ".vtu";
      string vtuFileName(ss.str());
      outFile << "        <DataSet index=\"" << i << "\" file=\"" << vtuFileName.c_str() << "\"></DataSet>" << std::endl;
    }
    outFile << "     </Block>" << std::endl;
  }
  outFile << "  </vtkMultiBlockDataSet>" << std::endl;
  outFile << "</VTKFile>" << std::endl;

  outFile.close();
}
