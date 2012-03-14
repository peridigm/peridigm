
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
  TEST_FOR_EXCEPTION( outputFormat != "BINARY",  std::invalid_argument, "PeridigmNS::OutputManager_ExodusII:::OutputManager_ExodusII() -- Output format must be BINARY for ExodusII.");

  // Default to not write full neighborlist
  writeNeighborlist = params->get<bool>("Bond Family",false); 
  TEST_FOR_EXCEPTION( (numProc != 1) && (writeNeighborlist),  std::invalid_argument, "PeridigmNS::OutputManager_ExodusII:::OutputManager_ExodusII() -- Parallel write of bond families not currently supported.");

  // Output filename base
  filenameBase = params->get<string>("Output Filename","dump"); 

  // User-requested fields for output 
  materialOutputFields = sublist(params, "Material Output Fields");

  // Initialize count (number of times write() has been called)
  // Initialize exodusCount (number of timesteps data actually written to exodus file)
  // Initialize to 0 because first call to write() corresponds to timestep 1
  exodusCount = count = 0;

  // Sentinal value for file handle
  file_handle = -1;

  // With ExodusII every object writes
  iWrite = true;

  // Rebalance count 
  rebalanceCount = 0;

  // Default to storing and writing doubles
  CPU_word_size = IO_word_size = sizeof(double);

  // Initialize the exodus database
  initializeExodusDatabase(blocks);

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
  exodusCount = exodusCount + 1;

  // Call compute manager; Updated any computed quantities before write
  peridigm->computeManager->compute(blocks);

  // Initialize new exodus database if needed
  // Each block is always rebalanced at the same time, so each datamanager should always return the same
  // rebalance count. Hence, we keep only a single static int for the rebalance count. If the first block 
  // rebalanced since last write, then all of them did. Force reinit database.
  if (rebalanceCount != blocks->begin()->getDataManager()->getRebalanceCount()) {
    rebalanceCount = blocks->begin()->getDataManager()->getRebalanceCount();
    initializeExodusDatabase(blocks);
  }

  // Open exodus database for writing
  float version;
  file_handle = ex_open(filename.str().c_str(), EX_WRITE, &CPU_word_size,&IO_word_size, &version);
  if (file_handle < 0) reportExodusError(file_handle, "write", "ex_open");

  // Write time value
  int retval = ex_put_time(file_handle,exodusCount,&current_time);
  if (retval!= 0) reportExodusError(retval, "write", "ex_put_time");
 
  int num_nodes = peridigm->getOneDimensionalMap()->NumMyElements();

  // Allocate temporary storage for all mothership-like data
  double *xptr = new double[num_nodes];
  double *yptr = new double[num_nodes];
  double *zptr = new double[num_nodes];

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
        double *block_ptr; block_ptr = NULL;
        // Exodus ignores element blocks when writing nodal variables
        if (fs.getRelation() == Field_ENUM::NODE) {
          // Loop over all blocks, copying data from each block into mothership-like vector
          std::vector<PeridigmNS::Block>::iterator blockIt;
          for(blockIt = blocks->begin(); blockIt != blocks->end() ; blockIt++) {
            Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
            Teuchos::RCP<Epetra_Vector> epetra_vector;
            if (fs.get_temporal() != Field_ENUM::TWO_STEP) // If stateless, get STEP_NONE
              epetra_vector = dataManager->getData(fs, Field_ENUM::STEP_NONE);
            else // If stateful, get STEP_NP1
              epetra_vector = dataManager->getData(fs, Field_ENUM::STEP_NP1);
            int block_num_nodes = (blockIt->getDataManager()->getOwnedScalarPointMap())->NumMyElements();
            epetra_vector->ExtractView(&block_ptr);
            // switch on dimension of data
            if (fs.getLength() == Field_ENUM::SCALAR) {
              // loop over contents of block vector; fill mothership-like vector
              for (int j=0;j<block_num_nodes; j++) {
                int GID = blockIt->getOwnedVectorPointMap()->GID(j);
                int msLID = peridigm->getOneDimensionalMap()->LID(GID);
                xptr[msLID] = block_ptr[j];
              }
            }
            else if (fs.getLength() == Field_ENUM::VECTOR2D) {
              // Peridgm doesn't have any 2D maps. 2D output not supported.
              TEST_FOR_EXCEPTION(true, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::write() -- 2D vector quantities not currently supported.");
            }
            else if (fs.getLength() == Field_ENUM::VECTOR3D) {
              // loop over contents of block vector; fill mothership-like vector
              for (int j=0;j<block_num_nodes; j++) {
                int GID = blockIt->getOwnedVectorPointMap()->GID(j);
                int msLID = peridigm->getThreeDimensionalMap()->LID(GID);
                xptr[msLID] = block_ptr[3*j];
                yptr[msLID] = block_ptr[3*j+1];
                zptr[msLID] = block_ptr[3*j+2];
              }
            } // end switch on data dimension
          } // end loop over blocks
          // Mothership-like vectors filled now; pass data to exodus database (switch again on dimension of data)
          if (fs.getLength() == Field_ENUM::SCALAR) {
            retval = ex_put_nodal_var(file_handle, exodusCount, node_output_field_map[nm], num_nodes, xptr);
            if (retval!= 0) reportExodusError(retval, "write", "ex_put_nodal_var");
          }
          else if (fs.getLength() == Field_ENUM::VECTOR2D) {
            // Peridgm doesn't have any 2D maps. 2D output not supported.
            TEST_FOR_EXCEPTION(true, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::write() -- 2D vector quantities not currently supported.");
          }
          else if (fs.getLength() == Field_ENUM::VECTOR3D) {
            // Writing all vector output as per-node data
            string tmpnameX = nm+"X";
            string tmpnameY = nm+"Y";
            string tmpnameZ = nm+"Z";
            retval = ex_put_nodal_var(file_handle, exodusCount, node_output_field_map[tmpnameX], num_nodes, xptr);
            if (retval!= 0) reportExodusError(retval, "write", "ex_put_nodal_var");
            retval = ex_put_nodal_var(file_handle, exodusCount, node_output_field_map[tmpnameY], num_nodes, yptr);
            if (retval!= 0) reportExodusError(retval, "write", "ex_put_nodal_var");
            retval = ex_put_nodal_var(file_handle, exodusCount, node_output_field_map[tmpnameZ], num_nodes, zptr);
            if (retval!= 0) reportExodusError(retval, "write", "ex_put_nodal_var");
          }
        } // end if per-node variable
        // Exodus wants element data written individually for each element block
        else if (fs.getRelation() == Field_ENUM::ELEMENT) {
          // Loop over all blocks, passing data from each block to exodus database
          std::vector<PeridigmNS::Block>::iterator blockIt;
          for(blockIt = blocks->begin(); blockIt != blocks->end() ; blockIt++) {
            int block_num_nodes = (blockIt->getDataManager()->getOwnedScalarPointMap())->NumMyElements();
            if (block_num_nodes == 0) continue; // Don't write data for empty blocks
            if (fs == Field_NS::GID) { // Handle special case of ID (int type)
              for (int j=0; j<block_num_nodes; j++)
                xptr[j] = (double)(((blockIt->getDataManager()->getOwnedScalarPointMap())->GID(j))+1);
              retval = ex_put_elem_var(file_handle, exodusCount, element_output_field_map[nm], blockIt->getID(), block_num_nodes, xptr);
              if (retval!= 0) reportExodusError(retval, "write", "ex_put_elem_var");
            }
            else if (fs == Field_NS::PROC_NUM) { // Handle special case of Proc_Num (int type)
              for (int j=0; j<block_num_nodes; j++)
                xptr[j] = (double)myPID;
              retval = ex_put_elem_var(file_handle, exodusCount, element_output_field_map[nm], blockIt->getID(), block_num_nodes, xptr);
              if (retval!= 0) reportExodusError(retval, "write", "ex_put_elem_var");
            }
            else {
              Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
              Teuchos::RCP<Epetra_Vector> epetra_vector;
              if (fs.get_temporal() != Field_ENUM::TWO_STEP) // If stateless, get STEP_NONE
                epetra_vector = dataManager->getData(fs, Field_ENUM::STEP_NONE);
              else // If stateful, get STEP_NP1
                epetra_vector = dataManager->getData(fs, Field_ENUM::STEP_NP1);
              epetra_vector->ExtractView(&block_ptr);
              // switch on dimension of data
              if (fs.getLength() == Field_ENUM::SCALAR) {
                retval = ex_put_elem_var(file_handle, exodusCount, element_output_field_map[nm], blockIt->getID(), block_num_nodes, block_ptr);
                if (retval!= 0) reportExodusError(retval, "write", "ex_put_elem_var");
              }
              else if (fs.getLength() == Field_ENUM::VECTOR2D) { 
                // Peridgm doesn't have any 2D maps. 2D output not supported.
                TEST_FOR_EXCEPTION(true, std::invalid_argument, "PeridigmNS::OutputManager_ExodusII::write() -- 2D vector quantities not currently supported.");
              }
              else if (fs.getLength() == Field_ENUM::VECTOR3D) {
                // separate out contents of block vector into x,y,z components
                for (int j=0;j<block_num_nodes; j++) {
                  int GID = blockIt->getOwnedVectorPointMap()->GID(j);
                  int msLID = peridigm->getThreeDimensionalMap()->LID(GID);
                  xptr[msLID] = block_ptr[3*j];
                  yptr[msLID] = block_ptr[3*j+1];
                  zptr[msLID] = block_ptr[3*j+2];
                }
                string tmpnameX = nm+"X";
                string tmpnameY = nm+"Y";
                string tmpnameZ = nm+"Z";
                retval = ex_put_nodal_var(file_handle, exodusCount, element_output_field_map[tmpnameX], block_num_nodes, xptr);
                if (retval!= 0) reportExodusError(retval, "write", "ex_put_nodal_var");
                retval = ex_put_nodal_var(file_handle, exodusCount, element_output_field_map[tmpnameY], block_num_nodes, yptr);
                if (retval!= 0) reportExodusError(retval, "write", "ex_put_nodal_var");
                retval = ex_put_nodal_var(file_handle, exodusCount, element_output_field_map[tmpnameZ], block_num_nodes, zptr);
                if (retval!= 0) reportExodusError(retval, "write", "ex_put_nodal_var");
              } // end switch on data dimension
            }
          } // end loop over blocks
        } // if per-element variable
      }
    }
  }

  delete[] xptr;
  delete[] yptr;
  delete[] zptr;

  // Flush write
  retval = ex_update(file_handle);
  if (retval!= 0) reportExodusError(retval, "write", "ex_update");
  retval = ex_close(file_handle);
  if (retval!= 0) reportExodusError(retval, "write", "ex_close");

}

// MLP: Change arguments to error reporting code

void PeridigmNS::OutputManager_ExodusII::initializeExodusDatabase(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks) {

  /*
   * First, determine name of output file
   */

  // Follow convention of replacing spaces or . with underscore
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
  // Construct output filename
  filename << filenameBase.c_str() << ".e";
  if (numProc > 1) {
    if (peridigm->analysisHasRebalance)
      filename << "-s" << setfill('0') << setw(5) << rebalanceCount;
    filename << "." << numProc;
    filename << "." << myPID;
  }

  /*
   * Now, initialize ExodusII database
   */

  // Default to storing and writing doubles
  int CPU_word_size, IO_word_size;
  CPU_word_size = IO_word_size = sizeof(double);

  // Initialize exodus database; Overwrite any existing file with this name
  file_handle = ex_create(filename.str().c_str(),EX_CLOBBER,&CPU_word_size,&IO_word_size);
  if (file_handle < 0) reportExodusError(file_handle, "OutputManager_ExodusII", "ex_create");

  // clear the maps
  element_output_field_map.clear();
  node_output_field_map.clear();

  // Initialize the database (assumes that Peridigm mothership vectors created and initialized)
  int num_dimensions = 3;
  int num_nodes = peridigm->getOneDimensionalMap()->NumMyElements();
  int num_elements = num_nodes;
  int num_element_blocks = blocks->size();
  int num_node_sets = 0;
  int num_side_sets = 0;
  // Initialize exodus file with parameters
  int retval = ex_put_init(file_handle,"Peridigm", num_dimensions, num_nodes, num_elements, num_element_blocks, num_node_sets, num_side_sets);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_init");

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
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_coord");

  // Write nodal coordinate names to database
  char *coord_names[3];
  coord_names[0] = "X";
  coord_names[1] = "Y";
  coord_names[2] = "Z";
  retval = ex_put_coord_names(file_handle,coord_names);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_coord_names");

  // Write element block parameters
  int *num_elem_in_block = new int[blocks->size()];
  int *num_nodes_in_elem = new int[blocks->size()];
  int *elem_block_ID     = new int[blocks->size()];
  std::vector<PeridigmNS::Block>::iterator blockIt;
  int i=0;
  for(i=0, blockIt = blocks->begin(); blockIt != blocks->end(); blockIt++, i++) {
    // Use only the number of owned elements
    num_elem_in_block[i] = (blockIt->getDataManager()->getOwnedScalarPointMap())->NumMyElements();
    num_nodes_in_elem[i] = 1; // always using sphere elements
    elem_block_ID[i]     = blockIt->getID();
    retval = ex_put_elem_block(file_handle,elem_block_ID[i],"SPHERE",num_elem_in_block[i],num_nodes_in_elem[i],0);
    if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_elem_block");
  }

  // Write element connectivity
  for(blockIt = blocks->begin(); blockIt != blocks->end(); blockIt++) {
    int numMyElements = blockIt->getOwnedScalarPointMap()->NumMyElements();
    if (numMyElements == 0) continue; // don't insert connectivity info for empty blocks
    int *connect = new int[numMyElements];
    for (int j=0;j<numMyElements;j++) {
      int GID = blockIt->getOwnedScalarPointMap()->GID(j);
      connect[j] = peridigm->getOneDimensionalMap()->LID(GID)+1;
    }
    retval = ex_put_elem_conn(file_handle, blockIt->getID(), connect);
    if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_elem_conn");
    delete[] connect;
  }

  // Write global node number map (global node IDs)
  int *node_map = new int[num_nodes];
  for (i=0; i<num_nodes; i++)
    node_map[i] = peridigm->getOneDimensionalMap()->GID(i)+1;
  retval = ex_put_node_num_map(file_handle, node_map);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_node_num_map");
  retval = ex_put_elem_num_map(file_handle, node_map);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_elem_num_map");

  // Create internal mapping of requested output fields to an integer.
  // The user requests output fields via strings, but Exodus wants an integer to index the output fields
  int node_output_field_index = 1;
  int element_output_field_index = 1;
  Teuchos::ParameterList::ConstIterator i1;
  // Loop over the material types in the materialOutputFields parameterlist
  for (i1 = materialOutputFields->begin(); i1 != materialOutputFields->end(); ++i1) {
    const Teuchos::ParameterEntry& val1 = materialOutputFields->entry(i1);
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
        if (fs.getLength() == Field_ENUM::SCALAR) {
          if (fs.getRelation() == Field_ENUM::NODE) {
            node_output_field_map.insert( std::pair<string,int>(nm,node_output_field_index) );
            node_output_field_index = node_output_field_index + 1;
          }
          else if (fs.getRelation() == Field_ENUM::ELEMENT) {
            element_output_field_map.insert( std::pair<string,int>(nm,element_output_field_index) );
            element_output_field_index = element_output_field_index + 1;
          }
        }
        if (fs.getLength() == Field_ENUM::VECTOR2D) {
          string tmpnameX = nm+"X";
          string tmpnameY = nm+"Y";
          if (fs.getRelation() == Field_ENUM::NODE) {
            node_output_field_map.insert( std::pair<string,int>(tmpnameX,node_output_field_index) );
            node_output_field_index = node_output_field_index + 1;
            node_output_field_map.insert( std::pair<string,int>(tmpnameY,node_output_field_index) );
            node_output_field_index = node_output_field_index + 1;
          }
          else if (fs.getRelation() == Field_ENUM::ELEMENT) {
            element_output_field_map.insert( std::pair<string,int>(tmpnameX,element_output_field_index) );
            element_output_field_index = element_output_field_index + 1;
            element_output_field_map.insert( std::pair<string,int>(tmpnameY,element_output_field_index) );
            element_output_field_index = element_output_field_index + 1;
          }
        }
        if (fs.getLength() == Field_ENUM::VECTOR3D) {
          string tmpnameX = nm+"X";
          string tmpnameY = nm+"Y";
          string tmpnameZ = nm+"Z";
          if (fs.getRelation() == Field_ENUM::NODE) {
            node_output_field_map.insert( std::pair<string,int>(tmpnameX,node_output_field_index) );
            node_output_field_index = node_output_field_index + 1;
            node_output_field_map.insert( std::pair<string,int>(tmpnameY,node_output_field_index) );
            node_output_field_index = node_output_field_index + 1;
            node_output_field_map.insert( std::pair<string,int>(tmpnameZ,node_output_field_index) );
            node_output_field_index = node_output_field_index + 1;
          }
          else if (fs.getRelation() == Field_ENUM::ELEMENT) {
            element_output_field_map.insert( std::pair<string,int>(tmpnameX,element_output_field_index) );
            element_output_field_index = element_output_field_index + 1;
            element_output_field_map.insert( std::pair<string,int>(tmpnameY,element_output_field_index) );
            element_output_field_index = element_output_field_index + 1;
            element_output_field_map.insert( std::pair<string,int>(tmpnameZ,element_output_field_index) );
            element_output_field_index = element_output_field_index + 1;
          }
        }
      }
    }
  }

  // Write information records

  // Write node var info 
  int num_node_vars = node_output_field_map.size();
  char **node_var_names = new char*[num_node_vars];
  for (i=0;i<num_node_vars;i++) node_var_names[i] = new char[MAX_STR_LENGTH+1]; // MAX_STR_LENGTH defined in ExodusII.h
  std::map<string,int>::iterator it;
  for ( it=node_output_field_map.begin() ; it != node_output_field_map.end(); it++ )
    strcpy(node_var_names[(it->second)-1], it->first.c_str() );
  retval = ex_put_var_param(file_handle,"N",num_node_vars);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_var_param");
  retval = ex_put_var_names(file_handle,"N",num_node_vars,node_var_names);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_var_names");

  // Write element var info 
  int num_element_vars = element_output_field_map.size();
  char **element_var_names = new char*[num_element_vars];
  for (i=0;i<num_element_vars;i++) element_var_names[i] = new char[MAX_STR_LENGTH+1]; // MAX_STR_LENGTH defined in ExodusII.h
  for ( it=element_output_field_map.begin() ; it != element_output_field_map.end(); it++ )
    strcpy(element_var_names[(it->second)-1], it->first.c_str() );
  retval = ex_put_var_param(file_handle,"E",num_element_vars);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_var_param");
  retval = ex_put_var_names(file_handle,"E",num_element_vars,element_var_names);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_var_names");

  // Write element truth table
  int *truthTable = new int[ blocks->size() * num_element_vars  ];
  for(i=0, blockIt = blocks->begin(); blockIt != blocks->end(); blockIt++)
    for (int j=0; j<num_element_vars; j++)
      truthTable[i++] = 1;
  retval = ex_put_elem_var_tab (file_handle, blocks->size(), num_element_vars, truthTable);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_put_var_tab");

  // Close file; re-open with call to write()
  retval = ex_update(file_handle);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_update");
  retval = ex_close(file_handle);
  if (retval!= 0) reportExodusError(retval, "initializeExodusDatabase", "ex_close");

  // Clean up
  delete[] xcoord_values; 
  delete[] ycoord_values; 
  delete[] zcoord_values;
  delete[] num_elem_in_block;
  delete[] num_nodes_in_elem;
  delete[] node_map;
  delete[] elem_block_ID;
  for (i = num_node_vars; i>0; i--) delete[] node_var_names[i-1];
  delete[] node_var_names;
  for (i = num_element_vars; i>0; i--) delete[] element_var_names[i-1];
  delete[] element_var_names;
}


void PeridigmNS::OutputManager_ExodusII::reportExodusError(int errorCode, const char *methodName, const char*exodusMethodName) {
  std::stringstream ss;
  if (errorCode < 0) { // error
    if (numProc > 1) ss << "Error on PID #" << myPID << ": ";
    ss << "PeridigmNS::OutputManager_ExodusII::" << methodName << "() -- Error code: " << errorCode << " (" << exodusMethodName << ")";
    TEST_FOR_EXCEPTION(1, std::invalid_argument, ss.str());
  }  
  else {
    if (numProc > 1) ss << "Warning on PID #" << myPID << ": ";
    ss << "PeridigmNS::OutputManager_ExodusII::" << methodName << "() -- Warning code: " << errorCode << " (" << exodusMethodName << ")";
    std::cout << ss.str() << std::endl;
  }
}


// MLP: Dead code that might be useful later

// Maybe use this code block later?
// It's more efficient to pull data from the mothership vectors, if they
// exist for the requested output data, rather than pulling from the
// vectors in the block
/*
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
*/


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

