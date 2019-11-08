/*! \file Peridigm_MeshConverter_Main.cpp
 *
 * Utility for converting hex/tet genesis mesh to a sphere-element mesh.
 */

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

#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
#endif

#include "Peridigm_HorizonManager.hpp"
#include "Peridigm_ExodusDiscretization.hpp"
#include "exodusII.h"
#include <iostream>
#include <sstream>
#include <Teuchos_YamlParameterListCoreHelpers.hpp>
#include <float.h>

void reportExodusError(int errorCode, const char *methodName, const char*exodusMethodName)
{
  std::stringstream ss;
  if (errorCode < 0) { // error
    ss << "MeshConverter -- Error code: " << errorCode << " (" << exodusMethodName << ")";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(1, ss.str());
  }
  else {
    ss << "MeshConverter -- Warning code: " << errorCode << " (" << exodusMethodName << ")";
    std::cout << ss.str() << std::endl;
  }
}

int main(int argc, char *argv[]) {

  // Banner
  std::cout << "\n-- Mesh Converter" << std::endl;

  // initialize MPI
  int mpi_id = 0;
  int mpi_size = 1;
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // comms
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  Teuchos::RCP<const Epetra_Comm> epetra_comm;
#ifdef HAVE_MPI
  epetra_comm = Teuchos::rcp(new Epetra_MpiComm(mpi_comm));
#else
  epetra_comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

  // input file
  if(argc != 2 && mpi_id == 0){
    std::cout << "Usage:  MeshConverter <input_file.yaml>\n" << std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  std::string yaml_file_name(argv[1]);
  Teuchos::RCP<Teuchos::ParameterList> params;
  params = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::Ptr<Teuchos::ParameterList> params_ptr(params.get());
  Teuchos::updateParametersFromYamlFile(yaml_file_name, params_ptr);

  Teuchos::ParameterList& mesh_conversion_params = params->sublist("Mesh Conversion", true);
  std::string output_file_name = mesh_conversion_params.get<std::string>("Output Mesh File");
  bool perform_neighborhood_search = mesh_conversion_params.get<bool>("Perform Neighborhood Search");
  std::string output_neighborhood_file_name = mesh_conversion_params.get<std::string>("Output Neighborhood File", "None");
  double scale_factor = mesh_conversion_params.get<double>("Scale Factor", 1.0);

  Teuchos::RCP<Teuchos::ParameterList> disc_params = Teuchos::rcpFromRef(params->sublist("Discretization", true));
  if(!disc_params->isParameter("Construct Neighborhood Lists")) {
    disc_params->set<bool>("Construct Neighborhood Lists", perform_neighborhood_search);
  }

  if (perform_neighborhood_search) {
    Teuchos::ParameterList& block_params = params->sublist("Blocks", true);
    PeridigmNS::HorizonManager::self().loadHorizonInformationFromBlockParameters(block_params);
  }

  PeridigmNS::ExodusDiscretization discretization(epetra_comm, disc_params);

  // Obtain the node sets
  Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets = discretization.getNodeSets();
  std::map< std::string, std::vector<int> >::iterator nsIt;

  int num_dimensions = 3;
  int num_elements = discretization.getNumElem();
  int num_nodes = num_elements;
  int num_element_blocks = discretization.getNumBlocks();
  int num_node_sets = nodeSets()->size();
  int num_side_sets = 0;

  // Add boundary layer nodes sets, if requested by user
  std::map< std::string, std::vector<int> > boundaryLayerNodeSets;
  if (mesh_conversion_params.isParameter("Boundary Layer Thickness")) {
    double boundaryLayerThickness = mesh_conversion_params.get<double>("Boundary Layer Thickness");

    // Determine the bounding box
    double x_min(DBL_MAX), x_max(-DBL_MAX), y_min(DBL_MAX), y_max(-DBL_MAX), z_min(DBL_MAX), z_max(-DBL_MAX);
    for( int i=0 ; i<num_elements ; i++ ) {
      int firstPoint = discretization.getInitialX()->Map().FirstPointInElement(i);
      double x = (*discretization.getInitialX())[firstPoint];
      double y = (*discretization.getInitialX())[firstPoint+1];
      double z = (*discretization.getInitialX())[firstPoint+2];
      if(x < x_min) x_min = x;
      if(x > x_max) x_max = x;
      if(y < y_min) y_min = y;
      if(y > y_max) y_max = y;
      if(z < z_min) z_min = z;
      if(z > z_max) z_max = z;
    }

    // Create boundary layer node sets
    boundaryLayerNodeSets["min_x_face"] = std::vector<int>();
    boundaryLayerNodeSets["max_x_face"] = std::vector<int>();
    boundaryLayerNodeSets["min_y_face"] = std::vector<int>();
    boundaryLayerNodeSets["max_y_face"] = std::vector<int>();
    boundaryLayerNodeSets["min_z_face"] = std::vector<int>();
    boundaryLayerNodeSets["max_z_face"] = std::vector<int>();
    boundaryLayerNodeSets["initial_velocity_node_set"] = std::vector<int>();
    for( int i=0 ; i<num_elements ; i++ ) {
      int firstPoint = discretization.getInitialX()->Map().FirstPointInElement(i);
      double x = (*discretization.getInitialX())[firstPoint];
      double y = (*discretization.getInitialX())[firstPoint+1];
      double z = (*discretization.getInitialX())[firstPoint+2];

      if(x < x_min + boundaryLayerThickness){
        boundaryLayerNodeSets["min_x_face"].push_back(i);
      }

      if(x > x_max - boundaryLayerThickness){
        boundaryLayerNodeSets["max_x_face"].push_back(i);
      }
      else {
        boundaryLayerNodeSets["initial_velocity_node_set"].push_back(i);
      }

      if(y < y_min + boundaryLayerThickness){
        boundaryLayerNodeSets["min_y_face"].push_back(i);
      }

      if(y > y_max - boundaryLayerThickness){
        boundaryLayerNodeSets["max_y_face"].push_back(i);
      }

      if(z < z_min + boundaryLayerThickness){
        boundaryLayerNodeSets["min_z_face"].push_back(i);
      }

      if(z > z_max - boundaryLayerThickness){
        boundaryLayerNodeSets["max_z_face"].push_back(i);
      }
    }
  }

  num_node_sets += boundaryLayerNodeSets.size();

  // Default to storing and writing doubles
  int CPU_word_size, IO_word_size;
  CPU_word_size = IO_word_size = sizeof(double);

  // Initialize exodus database; Overwrite any existing file with this name
  int file_handle = ex_create(output_file_name.c_str(), EX_CLOBBER, &CPU_word_size, &IO_word_size);
  if (file_handle < 0) reportExodusError(file_handle, "MeshConverter", "ex_create");

  // Initialize exodus file with parameters
  int retval = ex_put_init(file_handle, "MeshConverter", num_dimensions, num_nodes, num_elements, num_element_blocks, num_node_sets, num_side_sets);
  if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_init");
  //writeQARecord(file_handle);

  // Write the node sets
  unsigned int numNodeSets = nodeSets->size() + boundaryLayerNodeSets.size();
  int numNodesAcrossAllNodeSets = 0;
  for(nsIt = nodeSets->begin() ; nsIt != nodeSets->end() ; nsIt++)
    numNodesAcrossAllNodeSets += nsIt->second.size();
  for(nsIt = boundaryLayerNodeSets.begin() ; nsIt != boundaryLayerNodeSets.end() ; nsIt++)
    numNodesAcrossAllNodeSets += nsIt->second.size();
  std::vector<int> node_set_ids(numNodeSets);
  std::vector<int> num_nodes_per_set(numNodeSets);
  std::vector<int> num_dist_per_set(numNodeSets);
  std::vector<int> node_sets_node_index(numNodeSets);
  std::vector<int> node_sets_dist_index(numNodeSets);
  std::vector<int> node_sets_node_list(numNodesAcrossAllNodeSets);
  int* node_sets_dist_fact = 0;
  int nodeSetIndex = 0;
  int offset = 0;
  for(nsIt = nodeSets->begin() ; nsIt != nodeSets->end() ; nsIt++){
    std::vector<int>& nodeSet = nsIt->second;
    node_set_ids[nodeSetIndex] = nodeSetIndex + 1;
    num_nodes_per_set[nodeSetIndex] = nodeSet.size();
    num_dist_per_set[nodeSetIndex] = 0;
    node_sets_node_index[nodeSetIndex] = offset;
    node_sets_dist_index[nodeSetIndex] = 0;
    for(unsigned int i=0 ; i<nodeSet.size() ; ++i)
      node_sets_node_list[offset++] = nodeSet[i] + 1;
    nodeSetIndex += 1;
  }
  for(nsIt = boundaryLayerNodeSets.begin() ; nsIt != boundaryLayerNodeSets.end() ; nsIt++){
    std::vector<int>& nodeSet = nsIt->second;
    node_set_ids[nodeSetIndex] = nodeSetIndex + 1;
    num_nodes_per_set[nodeSetIndex] = nodeSet.size();
    num_dist_per_set[nodeSetIndex] = 0;
    node_sets_node_index[nodeSetIndex] = offset;
    node_sets_dist_index[nodeSetIndex] = 0;
    for(unsigned int i=0 ; i<nodeSet.size() ; ++i)
      node_sets_node_list[offset++] = nodeSet[i] + 1;
    nodeSetIndex += 1;
  }
  if(numNodeSets > 0){
    retval = ex_put_concat_node_sets(file_handle,
                                     &node_set_ids[0],
                                     &num_nodes_per_set[0],
                                     &num_dist_per_set[0],
                                     &node_sets_node_index[0],
                                     &node_sets_dist_index[0],
                                     &node_sets_node_list[0],
                                     node_sets_dist_fact);
    if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_concat_node_sets");
  }

  // Write the node set names
  char **node_set_names = NULL;
  if(numNodeSets > 0){
    node_set_names = new char*[numNodeSets];
    for(unsigned int i=0;i<numNodeSets;i++) node_set_names[i] = new char[MAX_STR_LENGTH+1];
    int index = 0;
    for(nsIt = nodeSets->begin() ; nsIt != nodeSets->end() ; nsIt++)
      strcpy(node_set_names[index++], nsIt->first.c_str());
    for(nsIt = boundaryLayerNodeSets.begin() ; nsIt != boundaryLayerNodeSets.end() ; nsIt++)
      strcpy(node_set_names[index++], nsIt->first.c_str());
    retval = ex_put_names(file_handle, EX_NODE_SET, node_set_names);
    if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_names EX_NODE_SET");
  }

  // Write nodal coordinate values
  // Exodus requires pointer to x,y,z coordinates of nodes, but Peridigm stores this data using a blockmap, which interleaves the data
  // So, extract and copy the data to temporary storage that can be handed to the exodus api
  double *coord_values;
  discretization.getInitialX()->ExtractView( &coord_values );
  int numMyElements = num_elements;
  std::vector<double> xcoord_values_vec(numMyElements), ycoord_values_vec(numMyElements), zcoord_values_vec(numMyElements);
  double *xcoord_values = &xcoord_values_vec[0];
  double *ycoord_values = &ycoord_values_vec[0];
  double *zcoord_values = &zcoord_values_vec[0];
  for( int i=0 ; i<numMyElements ; i++ ) {
    int firstPoint = discretization.getInitialX()->Map().FirstPointInElement(i);
    xcoord_values[i] = scale_factor*coord_values[firstPoint];
    ycoord_values[i] = scale_factor*coord_values[firstPoint+1];
    zcoord_values[i] = scale_factor*coord_values[firstPoint+2];
  }
  retval = ex_put_coord(file_handle,xcoord_values,ycoord_values,zcoord_values);
  if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_coord");

  // Write nodal coordinate names to database
  const char *coord_names[3] = {"x", "y", "z"};
  retval = ex_put_coord_names(file_handle,const_cast<char**>(coord_names));
  if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_coord_names");

  // Write element block parameters
  Teuchos::RCP< std::map< std::string, std::vector<int> > > blocks = discretization.getElementBlocks();
  std::vector<int> num_elem_in_block_vec(blocks->size()), num_nodes_in_elem_vec(blocks->size()), elem_block_ID_vec(blocks->size());
  int *num_elem_in_block = &num_elem_in_block_vec[0];
  int *num_nodes_in_elem = &num_nodes_in_elem_vec[0];
  int *elem_block_ID     = &elem_block_ID_vec[0];
  int num_attr = 2;
  std::map< std::string, std::vector<int> >::iterator blockIt;
  int i=0;
  for(i=0, blockIt = blocks->begin(); blockIt != blocks->end(); blockIt++, i++) {
    // Use only the number of owned elements
    num_elem_in_block[i] = blockIt->second.size();
    num_nodes_in_elem[i] = 1; // always using sphere elements
    elem_block_ID[i]     = discretization.blockNameToBlockId(blockIt->first);
    retval = ex_put_elem_block(file_handle, elem_block_ID[i], "SPHERE", num_elem_in_block[i], num_nodes_in_elem[i], num_attr);
    if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_elem_block");
  }

  // Write the block names
  TEUCHOS_TEST_FOR_EXCEPT_MSG(blocks->size() < 1, "\nMeshConverter, Zero element blocks found!\n");
  char **block_names = new char*[blocks->size()];
  for(unsigned int i=0 ; i<blocks->size() ; ++i)
    block_names[i] = new char[MAX_STR_LENGTH+1];
  int index = 0;
  for(blockIt = blocks->begin(); blockIt != blocks->end(); blockIt++)
    strcpy(block_names[index++], blockIt->first.c_str());
  retval = ex_put_names(file_handle, EX_ELEM_BLOCK, block_names);
  if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_names EX_ELEM_BLOCK");

  Teuchos::RCP<const Epetra_BlockMap> oneDimensionalMap = discretization.getGlobalOwnedMap(1);
  Teuchos::RCP<Epetra_Vector> cellVolume = discretization.getCellVolume();

  // Write element connectivity
  for(blockIt = blocks->begin(); blockIt != blocks->end(); blockIt++) {
    int numMyElements = blockIt->second.size();
    if (numMyElements == 0) continue; // don't insert connectivity info for empty blocks
    std::vector<int> connect_vec(numMyElements);
    int *connect = &connect_vec[0];
    std::vector<double> element_attributes(num_attr*numMyElements);
    for (int j=0;j<numMyElements;j++) {
      int globalId = blockIt->second[j];
      int localId = oneDimensionalMap->LID(globalId);
      connect[j] = localId + 1;
      double volume = (*cellVolume)[localId];
      double sphRadius = cbrt(volume);
      // by convention, the "sph radius" and nodal volume are stored as element attributes
      element_attributes[2*j] = sphRadius;
      element_attributes[2*j+1] = volume;
    }
    int blockId = discretization.blockNameToBlockId(blockIt->first);
    retval = ex_put_elem_conn(file_handle, blockId, connect);
    if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_elem_conn");
    retval = ex_put_elem_attr(file_handle, blockId, element_attributes.data());
    if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_elem_attr");
  }

  // Write global node number map (global node IDs)
  std::vector<int> node_map_vec(num_nodes);
  int *node_map = &node_map_vec[0];
  for (i=0; i<num_nodes; i++){
    node_map[i] = oneDimensionalMap->GID(i)+1;
  }
  retval = ex_put_node_num_map(file_handle, node_map);
  if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_node_num_map");

  // Write global element number map (global element IDs)
  std::vector<int> elem_map_vec(num_nodes);
  int *elem_map = &elem_map_vec[0];
  int elem_map_index = 0;
  for(blockIt = blocks->begin(); blockIt != blocks->end() ; blockIt++) {
    for(int i=0; i<blockIt->second.size() ; ++i){
      elem_map[elem_map_index++] = blockIt->second[i] + 1;
    }
  }
  retval = ex_put_elem_num_map(file_handle, elem_map);
  if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_put_elem_num_map");

  // Close file; re-open with call to write()
  retval = ex_update(file_handle);
  if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_update");
  retval = ex_close(file_handle);
  if (retval!= 0) reportExodusError(retval, "MeshConverter", "ex_close");

  if (perform_neighborhood_search) {
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhood_data = discretization.getNeighborhoodData();
    neighborhood_data->WriteToDisk(output_neighborhood_file_name);
  }

  // Clean up
  if(node_set_names != NULL){
    for (i = numNodeSets; i>0; i--) delete[] node_set_names[i-1];
    delete[] node_set_names;
  }
  if(block_names != NULL){
    for (unsigned int i = blocks->size(); i>0; i--) delete[] block_names[i-1];
    delete[] block_names;
  }

  if (mpi_id == 0) {
    std::cout << "\nSphere mesh written to:        " << output_file_name << std::endl;
    if (perform_neighborhood_search) {
      std::cout << "Neighborhood data written to:  " << output_neighborhood_file_name << std::endl;
    }
    std::cout << std::endl;
  }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  return 0;
}

