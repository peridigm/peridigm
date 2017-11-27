/*! \file Peridigm_InterfaceData.cpp */

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

#include <Peridigm_InterfaceData.hpp>
#include <exodusII.h>
#include <Epetra_Import.h>
#include <set>

void
PeridigmNS::InterfaceData::Initialize(std::vector<int> leftElements, std::vector<int> rightElements, std::vector<int> numNodesPerElem,
  std::vector<std::vector<int> > interfaceNodesVec, const Teuchos::RCP<const Epetra_Comm> & Comm){

  TEUCHOS_TEST_FOR_EXCEPTION(leftElements.size()!=rightElements.size(),std::invalid_argument,"");
  TEUCHOS_TEST_FOR_EXCEPTION(numNodesPerElem.size()!=rightElements.size(),std::invalid_argument,"");
  TEUCHOS_TEST_FOR_EXCEPTION(interfaceNodesVec.size()!=rightElements.size(),std::invalid_argument,"");
  comm = Comm;
  numOwnedPoints = leftElements.size();
  if(ownedIDs != 0)
    delete[] ownedIDs;
  ownedIDs = new int[numOwnedPoints];
  if(elementLeft != 0)
    delete[] elementLeft;
  if(elementRight != 0)
    delete[] elementRight;
  if(numNodes != 0)
    delete[] numNodes;
  elementLeft = new int[numOwnedPoints];
  elementRight = new int[numOwnedPoints];
  numNodes = new int[numOwnedPoints];
  interfaceMap = Teuchos::rcp(new Epetra_BlockMap(-1, numOwnedPoints, 1, 0, *comm));
  for(int i=0;i<numOwnedPoints;++i){
    const int interfaceGID = interfaceMap->GID(i);
    ownedIDs[i] = interfaceGID;
    elementLeft[i] = leftElements[i];
    elementRight[i] = rightElements[i];
    numNodes[i] = numNodesPerElem[i];
  }
  // initialize the fields on the interfaces
  interfaceAperture = Teuchos::rcp(new Epetra_Vector(*interfaceMap));
  interfaceAperture->PutScalar(0.0);

  // initialize storage of the nodes that make up each interface
  interfaceNodesMap = Teuchos::rcp(new Epetra_BlockMap(-1,numOwnedPoints,&ownedIDs[0],&numNodes[0],0,*comm));
  interfaceNodes = Teuchos::rcp(new Epetra_Vector(*interfaceNodesMap));
  interfaceNodes->PutScalar(-1.0);

  // copy over the nodes that make up each interface
  int elemIndex = 0;
  for(int i=0;i<numOwnedPoints;++i){
    elemIndex = interfaceNodes->Map().FirstPointInElement(i);
    for(int j=0;j<numNodesPerElem[i];++j){
      (*interfaceNodes)[elemIndex + j] = interfaceNodesVec[i][j];
    }
  }

  // generate an overlap map of the elements attached to all local interfaces:
  std::set<int> ElemGIDsSet;
  for(int i=0;i<numOwnedPoints;++i){
    ElemGIDsSet.insert(elementLeft[i]);
    ElemGIDsSet.insert(elementRight[i]);
  }
  const int numOverlapIds = ElemGIDsSet.size();

  int * overlapIds = new int[numOverlapIds];
  int index = 0;
  std::set<int>::iterator it;
  for(it=ElemGIDsSet.begin();it!=ElemGIDsSet.end();++it){
    overlapIds[index] = *it;
    index++;
  }

  elemOverlapMap = Teuchos::rcp(new Epetra_BlockMap(-1,numOverlapIds,&overlapIds[0],3,0,*comm));

  delete[] overlapIds;

}

void
PeridigmNS::InterfaceData::InitializeExodusOutput(Teuchos::RCP<Epetra_Vector> exodusMeshElementConnectivity, Teuchos::RCP<Epetra_Vector> exodusMeshNodePositions){

  if(comm->NumProc()>1){
    filename << "Interfaces.e." << comm->NumProc() << "." << comm->MyPID();
  }
  else{
    filename << "Interfaces.e";
  }
  std::string outputFileNameStr = filename.str();
  std::vector<char> writable(outputFileNameStr.size() + 1);
  std::copy(outputFileNameStr.begin(), outputFileNameStr.end(), writable.begin());

  int spaDim = 3; // force three dimensional output
  const int numShells = numOwnedPoints;
  const int numNodes = exodusMeshNodePositions->Map().NumMyElements();
  int error_int;
  int CPU_word_size = 0;
  int IO_word_size = 0;
  /* create EXODUS II file */
  const int output_exoid = ex_create (&writable[0],EX_CLOBBER,&CPU_word_size, &IO_word_size);
  exoid = output_exoid;

  // scan the connectivity to see if there are quads and tets:
  numQuads = 0;
  numTris = 0;
  for(int i=0;i<numOwnedPoints;++i){
    if(interfaceNodesMap->ElementSize(i)==4) numQuads++;
    if(interfaceNodesMap->ElementSize(i)==3) numTris++;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(numQuads+numTris!=numShells,std::logic_error,"numQuads " << numQuads << "  + numTris " << numTris << " should sum up to numShells " << numShells);


  const int numBlocks = 2;

  error_int = ex_put_init(exoid, &writable[0], spaDim, numNodes, numShells, numBlocks, 0, 0);
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"ex_put_init(): Failure");

  //  write initial coordinates and node/element maps
  float * x = new float[numNodes];
  float * y = new float[numNodes];
  float * z = new float[numNodes];

  int * shellMap = new int[numShells];
  int * nodeMap = new int[numNodes];

  for(int i=0;i<numNodes;++i){
    int nodeIndex = exodusMeshNodePositions->Map().FirstPointInElement(i);
    x[i] = (*exodusMeshNodePositions)[nodeIndex+0];
    y[i] = (*exodusMeshNodePositions)[nodeIndex+1];
    z[i] = (*exodusMeshNodePositions)[nodeIndex+2];
    nodeMap[i] = exodusMeshNodePositions->Map().GID(i) + 1; // numbering is one based in exodus
  }
  for(int i=0;i<numShells;++i){
    shellMap[i] = interfaceNodesMap->GID(i) +1; // numbering is one based in exodus
  }

  error_int = ex_put_coord(exoid, x, y, z);
  char * coord_names[3];
  coord_names[0] = (char*) "coordinates_x";
  coord_names[1] = (char*) "coordinates_y";
  coord_names[2] = (char*) "coordinates_x";
  error_int = ex_put_coord_names(exoid, coord_names);
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"ex_put_coord_names(): Failure");
  error_int = ex_put_elem_num_map(exoid, shellMap);
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"ex_put_elem_num_map(): Failure");
  error_int = ex_put_node_num_map(exoid, nodeMap);
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"ex_put_node_num_map(): Failure");
  delete[] shellMap;
  delete[] nodeMap;
  delete[] x; delete[] y; delete[] z;

  // write quad blocks:
  int blockIndex = 0;
  blockIndex ++;
  int num_nodes_per_elem_q4 = 4;
  const std::string elem_type_str_q4 = numQuads > 0 ? "QUAD4" : "NULL";
  char * elem_type_q4 = const_cast<char *>(elem_type_str_q4.c_str());
  const int numElemInBlock_q4 = numQuads;
  error_int = ex_put_elem_block(exoid, blockIndex, elem_type_q4, numElemInBlock_q4, num_nodes_per_elem_q4, 0); // no attributes put in output file
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"ex_put_elem_block(): Failure");
  blockIndex ++;
  int num_nodes_per_elem_t3 = 3;
  const std::string elem_type_str_t3 = numTris > 0 ? "TRI3": "NULL";
  char * elem_type_t3 = const_cast<char *>(elem_type_str_t3.c_str());
  const int numElemInBlock_t3 = numTris;
  error_int = ex_put_elem_block(exoid, blockIndex, elem_type_t3, numElemInBlock_t3, num_nodes_per_elem_t3, 0); // no attributes put in output file
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"ex_put_elem_block(): Failure");

  //  write elem connectivities
  //  HEADS UP: the connectivities will not write to file until ex_close is called
  //  also note that the nodes must be the local node ids

  const int connLength = interfaceNodes->Map().NumMyElements();

  blockIndex = 0;
  blockIndex ++;
  int conn_index = 0;
  int * block_connect_q4 = new int[numQuads*num_nodes_per_elem_q4];
  for(int it=0;it<connLength;++it)
  {
    if(interfaceNodes->Map().ElementSize(it)!=4) continue;
    int elemIndex = interfaceNodes->Map().FirstPointInElement(it);
    for(int nn=0;nn<4;++nn){
      int node = static_cast<int>( (*interfaceNodes)[elemIndex+nn] );
      block_connect_q4[conn_index*num_nodes_per_elem_q4+nn] = exodusMeshNodePositions->Map().LID(node) + 1;// nodes are 1 based in exodus
    }
    conn_index++;
  }
  error_int = ex_put_elem_conn(exoid, blockIndex, block_connect_q4);
  delete[] block_connect_q4;
  blockIndex ++;
  conn_index = 0;
  int * block_connect_t3 = new int[numTris*num_nodes_per_elem_t3];
  for(int it=0;it<connLength;++it)
  {
    if(interfaceNodes->Map().ElementSize(it)!=3) continue;
    int elemIndex = interfaceNodes->Map().FirstPointInElement(it);
    for(int nn=0;nn<3;++nn){
      int node = static_cast<int>( (*interfaceNodes)[elemIndex+nn] );
      block_connect_t3[conn_index*num_nodes_per_elem_t3+nn] = exodusMeshNodePositions->Map().LID(node) + 1;// nodes are 1 based in exodus
    }
    conn_index++;
  }
  error_int = ex_put_elem_conn(exoid, blockIndex, block_connect_t3);
  delete[] block_connect_t3;

  int numVariables = 1; // TODO careful with this, if more fields are added to the interface data this must be updated
  char** eleVarNames = new char*[numVariables];
  std::vector<std::string> strNames;
  strNames.push_back("interface_aperture");
  for (int i=0;i<numVariables;++i)
    eleVarNames[i] = (char*) (strNames[i].c_str());
  error_int = ex_put_var_param(exoid, (char*) "e", numVariables);
  error_int = ex_put_var_names(exoid, (char*) "e", numVariables, &eleVarNames[0]);

  delete [] eleVarNames;

  // write the truth table
  int * truth_tab = new int[numVariables*numBlocks];
  for(int i=0;i<numVariables*numBlocks;++i)
    truth_tab[i] = 1;
  error_int = ex_put_elem_var_tab (exoid, numBlocks, numVariables, truth_tab);
  delete [] truth_tab;

  error_int = ex_update(exoid);
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"Exodus file close failed.");
  error_int = ex_close(exoid);
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"Exodus file close failed.");

}

void
PeridigmNS::InterfaceData::WriteExodusOutput(int timeStep, const float & timeValue, Teuchos::RCP<Epetra_Vector> x, Teuchos::RCP<Epetra_Vector> y){

  int error_int = 0;

  int CPU_word_size = 0;
  int IO_word_size = 0;
  float version = 0;
  std::string outputFileNameStr = filename.str();
  std::vector<char> writable(outputFileNameStr.size() + 1);
  std::copy(outputFileNameStr.begin(), outputFileNameStr.end(), writable.begin());

  exoid = ex_open(&writable[0], EX_WRITE, &CPU_word_size, &IO_word_size, &version);

  error_int = ex_put_time(exoid, timeStep, &timeValue);
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error, "ex_put_time(): Failure");

  float * quadValues = new float[numQuads];
  float * triValues = new float[numTris];

  // populate the quad values
  int quadIndex = 0;
  int triIndex = 0;
  for(int i=0;i<numOwnedPoints;++i){
    if(interfaceNodesMap->ElementSize(i)==4){
      quadValues[quadIndex] = (*interfaceAperture)[i];
      quadIndex++;
    }
    else if(interfaceNodesMap->ElementSize(i)==3){
      triValues[triIndex] = (*interfaceAperture)[i];
      triIndex++;
    }
    else{
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"size of this element is not recognized: " << interfaceNodesMap->ElementSize(i));
    }
  }

  int blockIndex = 0;
  const int varIndex = 1;
  blockIndex++;
  if(numQuads > 0){
    error_int = ex_put_elem_var(exoid, timeStep, varIndex, blockIndex, numQuads, &quadValues[0]);
    TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"Failure ex_put_elem_var(): ");
  }
  blockIndex++;
  if(numTris > 0){
    error_int = ex_put_elem_var(exoid, timeStep, varIndex, blockIndex, numTris, &triValues[0]);
    TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"Failure ex_put_elem_var(): ");
  }

  delete [] quadValues;
  delete [] triValues;

  // update the apertures...
  // import the mothership vectors x and y to the overlap epetra vectors
  Teuchos::RCP<const Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*elemOverlapMap, x->Map()));

  Teuchos::RCP<Epetra_Vector> xOverlap = Teuchos::rcp(new Epetra_Vector(*elemOverlapMap,true));
  xOverlap->Import(*x,*importer,Insert);
  Teuchos::RCP<Epetra_Vector> yOverlap = Teuchos::rcp(new Epetra_Vector(*elemOverlapMap,true));
  yOverlap->Import(*y,*importer,Insert);

  double *xValues;
  xOverlap->ExtractView( &xValues );
  double *yValues;
  yOverlap->ExtractView( &yValues );

  double xLeft=0,yLeft=0,zLeft=0,xRight=0,yRight=0,zRight=0;
  double XLeft=0,YLeft=0,ZLeft=0,XRight=0,YRight=0,ZRight=0;
  double X=0,Y=0;
  double dx=0,dy=0,dz=0,dX=0,dY=0,dZ=0;
  int elemIndexLeft=-1,elemIndexRight=-1,GIDLeft=-1,GIDRight=-1;

  for(int i=0;i<numOwnedPoints;++i){
    GIDLeft = elementLeft[i];
    GIDRight = elementRight[i];

    elemIndexLeft = xOverlap->Map().FirstPointInElement(elemOverlapMap->LID(GIDLeft));
    elemIndexRight = xOverlap->Map().FirstPointInElement(elemOverlapMap->LID(GIDRight));

    xLeft = xValues[elemIndexLeft+0];
    yLeft = xValues[elemIndexLeft+1];
    zLeft = xValues[elemIndexLeft+2];
    xRight = xValues[elemIndexRight+0];
    yRight = xValues[elemIndexRight+1];
    zRight = xValues[elemIndexRight+2];

    XLeft = yValues[elemIndexLeft+0];
    YLeft = yValues[elemIndexLeft+1];
    ZLeft = yValues[elemIndexLeft+2];
    XRight = yValues[elemIndexRight+0];
    YRight = yValues[elemIndexRight+1];
    ZRight = yValues[elemIndexRight+2];

    dx = xRight - xLeft;
    dy = yRight - yLeft;
    dz = zRight - zLeft;

    dX = XRight - XLeft;
    dY = YRight - YLeft;
    dZ = ZRight - ZLeft;

    X = std::sqrt(dx*dx + dy*dy + dz*dz);
    Y = std::sqrt(dX*dX + dY*dY + dZ*dZ);

    interfaceAperture->ReplaceMyValue(i,0,Y-X);
  }
  error_int = ex_update(exoid);
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"Exodus file close failed.");
  error_int = ex_close(exoid);
  TEUCHOS_TEST_FOR_EXCEPTION(error_int,std::logic_error,"Exodus file close failed.");

}
