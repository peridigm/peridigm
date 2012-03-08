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

#include <netcdf.h>
#include <exodusII.h>

#include <Epetra_Comm.h>
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include <Teuchos_TestForException.hpp>

#include "Peridigm.hpp"
#include "Peridigm_OutputManager_ExodusII.hpp"
#include "mesh_output/Field.h"


PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII(const Teuchos::RCP<Teuchos::ParameterList>& params, 
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
  if (!isValid) TEST_FOR_EXCEPTION(1, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII:::OutputManager_ExodusII() -- Invalid parameter, type or value.");

  try {
    numProc = params->INVALID_TEMPLATE_QUALIFIER get<int>("NumProc");
  }
  catch ( const std::exception& e) {
    TEST_FOR_EXCEPTION(1, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII:::OutputManager_ExodusII() -- numProc not present.");
  }

  try {
    myPID = params->INVALID_TEMPLATE_QUALIFIER get<int>("MyPID");
  }
  catch ( const std::exception& e) {
    TEST_FOR_EXCEPTION(1,  std::invalid_argument, "PeridigmNS::OutputManager_ExodusII:::OutputManager_ExodusII() -- MyPID not present.");
  }

  // Default to no output
  frequency = params->get<int>("Output Frequency",-1); 

  // Default to BINARY output
  outputFormat = params->get<string>("Output Format","BINARY"); 
  TEST_FOR_EXCEPTION( (outputFormat != "ASCII") && (outputFormat != "BINARY"),  std::invalid_argument, "PeridigmNS::OutputManager_ExodusII:::OutputManager_ExodusII() -- Unknown output format. Must be ASCII or BINARY.");

  // Default to not write full neighborlist
  writeNeighborlist = params->get<bool>("Bond Family",false); 
  TEST_FOR_EXCEPTION( (numProc != 1) && (writeNeighborlist),  std::invalid_argument, "PeridigmNS::OutputManager_ExodusII:::OutputManager_ExodusII() -- Parallel write of bond families not currently supported.");

  // Output filename base
  filenameBase = params->get<string>("Output Filename","dump"); 

  // User-requested fields for output 
  materialOutputFields = sublist(params, "Material Output Fields");

  // Initialize count (number of times write() has been called)
  // Initialize exodus_count (number of timesteps data actually written to exodus file)
  // Initialize to 0 because first call to write() corresponds to timestep 1
  exodus_count = count = 0;

  // Sentinal value for file handle
  file_handle = -1;

  // With ExodusII every object writes
  iWrite = true;

  // Create container to hold processor ID # data
  proc_num = Teuchos::rcp( new std::vector<int>() );

  // VTK doesn't like spaces or . so replace them with underscore
  {
  int warningFlag = 0;
  string outString;
  outString.append("\n\n***WARNING***\n");
  outString.append("PeridigmNS::OutputManager_ExodusII:::OutputManager_ExodusII() -- Avoid use of filenames containing '.' (period) and ' ' (space) with ExodusII.\n"); 
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

  // MLP: Todo - Check that response to return codes is valid in all cases

  // Construct output filename
  std::ostringstream filename; 
  filename << filenameBase.c_str() << ".e";
  if (numProc > 1) {
    filename << "." << numProc;
    filename << "." << myPID;
  }
  filename.str().c_str();

  // Default to storing and writing doubles
  int CPU_word_size, IO_word_size;
  CPU_word_size = IO_word_size = sizeof(double);

  // Initialize exodus database; Overwrite any existing file with this name
  file_handle = ex_create(filename.str().c_str(),EX_CLOBBER,&CPU_word_size,&IO_word_size);
  TEST_FOR_EXCEPTION(file_handle<0, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Invalid output filename.");

  // Initialize the database (assumes that Peridigm mothership vectors created and initialized)
  int num_dimensions = 3;
  int num_nodes = peridigm->x->Map().NumMyElements();
  int num_elements = num_nodes;
  int num_element_blocks = blocks->size();
  int num_node_sets = 0;
  int num_side_sets = 0;
  // Initialize exodus file with parameters
  int retval = ex_put_init(file_handle,"Peridigm", num_dimensions, num_nodes, num_elements, num_element_blocks, num_node_sets, num_side_sets);
  TEST_FOR_EXCEPTION(retval<0, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Error opening output ExodusII database (ex_put_init()).");
  if (retval > 0) std::cout << "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Warning upon opening output ExodusII database (ex_put_init()). Warning code: " << retval << std::endl;

  // Write nodal coordinate values
  // Exodus requires pointer to x,y,z coordinates of nodes, but Peridigm stores this data using a blockmap, which interleaves the data 
  // So, extract and copy the data to temporary storage that can be handed to the exodus api
  double *coord_values;
  peridigm->x->ExtractView( &coord_values );
  int numMyElements = peridigm->x->Map().NumMyElements();
  double *xcoord_values = new double[numMyElements];
  double *ycoord_values = new double[numMyElements];
  double *zcoord_values = new double[numMyElements];
  for( int i=0 ; i<numMyElements ; i++ ) {
    int firstPoint = peridigm->x->Map().FirstPointInElement(i);
    xcoord_values[i] = coord_values[firstPoint];
    ycoord_values[i] = coord_values[firstPoint+1];
    zcoord_values[i] = coord_values[firstPoint+2];
  }
  retval = ex_put_coord(file_handle,xcoord_values,ycoord_values,zcoord_values);
  TEST_FOR_EXCEPTION(retval<0, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Error opening output ExodusII database (ex_put_init()).");
  if (retval > 0) std::cout << "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Warning upon opening output ExodusII database (ex_put_init()). Warning code: " << retval << std::endl;

  // Write nodal coordinate names to database
  char *coord_names[3];
  coord_names[0] = "X";
  coord_names[1] = "Y";
  coord_names[2] = "Z";
  retval = ex_put_coord_names(file_handle,coord_names);
  TEST_FOR_EXCEPTION(retval<0, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Error opening output ExodusII database.");
  if (retval > 0) std::cout << "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Warning upon opening output ExodusII database. Warning code: " << retval << std::endl;

  // Write element block parameters
  int *num_elem_in_block = new int[blocks->size()];
  int *num_nodes_in_elem = new int[blocks->size()];
  int *elem_block_ID     = new int[blocks->size()];
  std::vector<PeridigmNS::Block>::iterator blockIt;
  int i=0;
  for(i=0, blockIt = blocks->begin(); blockIt != blocks->end(); blockIt++, i++) {
    double *xptr;
    Teuchos::RCP<Epetra_Vector> myX =  blockIt->getDataManager()->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE);
    myX->ExtractView( &xptr );
    // Use only the number of owned elements
    num_elem_in_block[i] = (blockIt->getDataManager()->getOwnedScalarPointMap())->NumMyElements();
    num_nodes_in_elem[i] = 1; // always using sphere elements
    elem_block_ID[i]     = blockIt->getID();
    retval = ex_put_elem_block(file_handle,elem_block_ID[i],"SPHERE",num_elem_in_block[i],num_nodes_in_elem[i],0);
    TEST_FOR_EXCEPTION(retval<0, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Error opening output ExodusII database.");
    if (retval > 0) std::cout << "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Warning upon opening output ExodusII database. Warning code: " << retval << std::endl;
  }

  // Write element connectivity
  for(i=0, blockIt = blocks->begin(); blockIt != blocks->end(); blockIt++, i++) {
    int numMyElements = blockIt->getOwnedScalarPointMap()->NumMyElements();
    int *connect = new int[numMyElements];
    for (int j=0;j<numMyElements;j++) {
      int GID = blockIt->getOwnedScalarPointMap()->GID(j);
      connect[j] = peridigm->getOneDimensionalMap()->LID(GID)+1;
    }
    retval = ex_put_elem_conn(file_handle, elem_block_ID[i], connect);
    TEST_FOR_EXCEPTION(retval<0, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Error opening output ExodusII database.");
    if (retval > 0) std::cout << "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Warning upon opening output ExodusII database. Warning code: " << retval << std::endl;
    delete[] connect;
  }

  // Write global node number map (global node IDs)
  int *node_map = new int[num_nodes];
  for (i=0; i<num_nodes; i++)
    node_map[i] = peridigm->getOneDimensionalMap()->GID(i)+1;
  retval = ex_put_node_num_map(file_handle, node_map);
  TEST_FOR_EXCEPTION(retval<0, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Error opening output ExodusII database.");
  if (retval > 0) std::cout << "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Warning upon opening output ExodusII database. Warning code: " << retval << std::endl;
  retval = ex_put_elem_num_map(file_handle, node_map);
  TEST_FOR_EXCEPTION(retval<0, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Error opening output ExodusII database.");
  if (retval > 0) std::cout << "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Warning upon opening output ExodusII database. Warning code: " << retval << std::endl;

  // Write information records
  // MLP: For now, just hardcode to displacement (x,y,z)
  int num_nodal_vars = 3;
  char *var_names[3];
  var_names[0] = "DisplacementX";
  var_names[1] = "DisplacementY";
  var_names[2] = "DisplacementZ";
  retval = ex_put_var_param(file_handle,"N",num_nodal_vars);
  retval = ex_put_var_names(file_handle,"N",num_nodal_vars,var_names);

  // Close file; re-open with call to write()
  retval = ex_update(file_handle);
//  retval = ex_close(file_handle);

  // Clean up
  delete[] xcoord_values; 
  delete[] ycoord_values; 
  delete[] zcoord_values;
  delete[] num_elem_in_block;
  delete[] num_nodes_in_elem;
  delete[] node_map;
  delete[] elem_block_ID;

}

Teuchos::ParameterList PeridigmNS::OutputManager_ExodusII::getValidParameterList() {

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

PeridigmNS::OutputManager_ExodusII::~OutputManager_ExodusII() {

}

void PeridigmNS::OutputManager_ExodusII::write(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks, double current_time) {

  if (!iWrite) return;

  // increment index count
  count = count + 1;

  // Only write if frequency count match
  if (frequency<=0 || count%frequency!=0) return;

  // increment exodus_count index
  exodus_count = exodus_count + 1;

  // Call compute manager; Updated any computed quantities before write
  peridigm->computeManager->compute(blocks);

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

  // MLP: What to do with Exodus database after rebalance? Exodus assumes that 
  // the model is described by data which are static (do not change through time). This data
  // includes nodal coordinates, element connectivity (node lists for each element), element
  // attributes, and node sets and side sets (used to aid in applying loading conditions and
  // boundary constraints).

/*
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
*/

//  if (blocks->size() == 1) // if only one block, don't write with block IDs
//    this->write(blocks->begin()->getDataManager(),blocks->begin()->getNeighborhoodData(),grids.front(),vtkWriters.front(),current_time);
//  else { // there are many blocks; write output using block IDs
//    std::vector< vtkSmartPointer<vtkExodusIIWriter>  >::iterator writerIt;
//    for(blockIt = blocks->begin(), gridIt=grids.begin(), writerIt=vtkWriters.begin() ; blockIt != blocks->end() ; blockIt++, gridIt++, writerIt++)
//      this->write(blockIt->getDataManager(),blockIt->getNeighborhoodData(),*gridIt,*writerIt,current_time,blockIt->getID());
//  }

  // Write time value
  int retval = ex_put_time(file_handle,exodus_count,&current_time);
  TEST_FOR_EXCEPTION(retval<0, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Error opening output ExodusII database (ex_put_init()).");
  if (retval > 0) std::cout << "PeridigmNS::OutputManager_ExodusII::OutputManager_ExodusII() -- Warning upon opening output ExodusII database (ex_put_init()). Warning code: " << retval << std::endl;
 
  // MLP: Only displacement for block #1 right now.
  //int ex_put_nodal_var (exoid, time_step, nodal_var_index,num_nodes, nodal_var_vals);
  int num_nodes = peridigm->u->Map().NumMyElements();

  int numMyElements = peridigm->x->Map().NumMyElements();
  double *xdisp = new double[numMyElements];
  double *ydisp = new double[numMyElements];
  double *zdisp = new double[numMyElements];
  double *disp_values;
  peridigm->u->ExtractView( &disp_values );
  for( int i=0 ; i<numMyElements ; i++ ) {
    //int firstPoint = peridigm->u->Map().FirstPointInElement(i);
    int firstPoint = 3*i;
    xdisp[i] = disp_values[firstPoint];
    ydisp[i] = disp_values[firstPoint+1];
    zdisp[i] = disp_values[firstPoint+2];
  }
  retval = ex_put_nodal_var(file_handle, exodus_count, 1, num_nodes, xdisp);
  retval = ex_put_nodal_var(file_handle, exodus_count, 2, num_nodes, ydisp);
  retval = ex_put_nodal_var(file_handle, exodus_count, 3, num_nodes, zdisp);
   

/*
  for(blockIt = blocks->begin(), gridIt=grids.begin() ; blockIt != blocks->end() ; blockIt++, gridIt++) {
     // step #1: associate requested field data with each grid object

     //this->write(blockIt->getDataManager(),blockIt->getNeighborhoodData(),*gridIt,*writerIt,current_time,blockIt->getID());

     Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
     Teuchos::RCP<const NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
     vtkSmartPointer<vtkUnstructuredGrid> grid = *gridIt;

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

     // step #2: hand pointer to grid object to multi block data set
     mbds->SetBlock(blockIt->getID(), grid);

  }
*/

  // *********************
  // setup model meta data
  // *********************
/*
  metadata = vtkSmartPointer<vtkModelMetadata>::New();
  metadata->SetTitle("Title");
  metadata->SetNumberOfNodeSets(0); 
  metadata->SetNumberOfSideSets(0); 
  metadata->SetNumberOfBlocks(blocks->size());

  char **names = new char*[3];
  for(int i=0;i<3;i++) { names[i] = new char[2]; };
  names[0][0] = 'X'; names[0][1] = '\0';
  names[1][0] = 'Y'; names[1][1] = '\0';
  names[2][0] = 'Z'; names[2][1] = '\0';
  metadata->SetCoordinateNames  (3, names);

  // see vtkExodusIIWriter CreateBlockIDMetadata
  int *blockIds = new int[blocks->size()];
  char **blockNames = new char *[blocks->size()];
  int *numElements = new int[blocks->size()];
  int *numNodesPerElement = new int[blocks->size()];
  int *numAttributes = new int[blocks->size()];

  int i=0;
  for(i=0, blockIt = blocks->begin(); blockIt != blocks->end() ; i++, blockIt++) { 
    // Set block ID
    blockIds[i] = blockIt->getID(); 
    // Set block name
    char *nm = new char[12];
    strcpy(nm, "sphere"); // see vtkExodusIIWriter::GetCellTypeName(), VTK_VERTEX maps to "sphere"
    blockNames[i] = nm;
    // Set number of elements in block
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<Epetra_Vector> myX = dataManager->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE);
    const Epetra_BlockMap& xMap = myX->Map();
    numElements[i] = xMap.NumMyElements();
    // Set number of nodes per element (always 1)
    numNodesPerElement[i] = 1;
    // Set number of attributes per element
    numAttributes[i] = 0;
  }
  metadata->SetBlockIds(blockIds);
  metadata->SetBlockElementType(blockNames);
  metadata->SetBlockNumberOfElements(numElements);
  metadata->SetBlockNodesPerElement(numNodesPerElement);
  metadata->SetBlockNumberOfAttributesPerElement(numAttributes);

  // Set global variable names (none)
  char **globalVariableNames = NULL;
  metadata->SetGlobalVariableNames(0, globalVariableNames);

  float *times = new float[11];
  times[0] = 0.0; times[1] = 0.1; times[2] = 0.2; times[3] = 0.3;
  times[4] = 0.4; times[5] = 0.5; times[6] = 0.6; times[7] = 0.7;
  times[8] = 0.8; times[9] = 0.9; times[10] = 1.0;
  metadata->SetTimeSteps(11, times);

  static int timestep = 0;
  timestep++;
  metadata->SetTimeStepIndex(timestep);

//  Need to embed cell array in unstructured grid data?
//  metadata->SetBlockIdArrayName  ( const char *   ) [virtual] 

  vtkWriter->SetModelMetadata(metadata);

  // metadata->PrintGlobalInformation();

  // you must SetBlockIds and SetBlockIdArrayName. Your vtkUnstructuredGrid must have a cell array giving the block ID for each cell.

  // SetInput for writer as MBDS, then call write
  //vtkWriter->SetInput(mbds);

  vtkSmartPointer<vtkTrivialProducer> tp = vtkSmartPointer<vtkTrivialProducer>::New();
  tp->SetOutput(mbds);
  vtkWriter->SetInputConnection(tp->GetOutputPort());

  vtkWriter->Update();
  vtkWriter->Write();

exit(1);
*/

  // Flush write
  retval = ex_update(file_handle);

}

/*
void PeridigmNS::OutputManager_ExodusII::write(Teuchos::RCP<PeridigmNS::DataManager> dataManager,
                                              Teuchos::RCP<const NeighborhoodData> neighborhoodData,
                                              vtkSmartPointer<vtkUnstructuredGrid> grid,
                                              vtkSmartPointer<vtkExodusIIWriter> vtkWriter,
                                              double current_time, int block_id) {

  return;

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
  //vtkWriter->writeTimeStep(current_time,grid,block_id);
  vtkWriter->Write();
}
*/

