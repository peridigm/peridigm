/*! \file Peridigm_ZoltanSearchTree.cpp */
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

#include "Peridigm_ZoltanSearchTree.hpp"
// #include <stdexcept>
// #include <stdlib.h>
// #include <stdio.h>
#include <Teuchos_Assert.hpp>


PeridigmNS::ZoltanSearchTree::ZoltanSearchTree(int numPoints, double* coordinates)
  : SearchTree(numPoints, coordinates),
    callbackdata(numPoints,coordinates),zoltan(0),
    partToCoordIdx(new int[numPoints]),
    searchParts(new int[numPoints])
{



  zoltan = Zoltan_Create(MPI_COMM_SELF);
  Zoltan_Set_Param(zoltan, "debug_level", "0");
  Zoltan_Set_Param(zoltan, "rcb_output_level", "0");
  /*
   * query function returns the number of objects that are currently assigned to the processor
   */
  Zoltan_Set_Num_Obj_Fn(zoltan, &get_num_points, &callbackdata);

  /*
   * query function fills two (three if weights are used) arrays with information about the objects
   * currently assigned to the processor. Both arrays are allocated (and subsequently freed) by Zoltan;
   * their size is determined by a call to a ZOLTAN_NUM_OBJ_FN query function to get the array size.
   */
  Zoltan_Set_Obj_List_Fn(zoltan, &get_point_ids, &callbackdata);

  /*
   * query function returns the number of values needed to express the geometry of an object.
   * For example, for a two-dimensional mesh-based application, (x,y) coordinates are needed
   * to describe an object's geometry; thus the ZOLTAN_NUM_GEOM_FN query function should return
   * the value of two. For a similar three-dimensional application, the return value should be three.
   */
  Zoltan_Set_Num_Geom_Fn(zoltan, &get_dimension, NULL);

  /*
   * query function returns a vector of geometry values for a list of given objects. The geometry
   * vector is allocated by Zoltan to be of size num_obj * num_dim.
   */

  Zoltan_Set_Geom_Multi_Fn(zoltan, &get_point_coordinates, &callbackdata);

  /*
   * setting method
   */
  Zoltan_Set_Param(zoltan,"lb_method","rcb");

  /*
   * keep cuts
   */
  Zoltan_Set_Param(zoltan,"keep_cuts","1");

  Zoltan_Set_Param(zoltan,"return_lists","part");
  Zoltan_Set_Param(zoltan, "num_lid_entries", "0");
  char s[11];
  sprintf(s,"%d",numPoints);
  Zoltan_Set_Param(zoltan, "num_global_parts", s);

  int numimp = -1, numexp = -1;
  int numgid = 1, numlid = 0;
  ZOLTAN_ID_PTR impgid=NULL, implid=NULL;
  ZOLTAN_ID_PTR gid=NULL, lid=NULL;
  int *imppart = NULL, *impproc = NULL;
  int *procs = NULL, *parts=NULL;


  int changes;

  int ierr = Zoltan_LB_Partition(zoltan, &changes, &numgid, &numlid,
                                 &numimp, &impgid, &implid, &imppart, &impproc,
                                 &numexp, &gid, &lid, &procs, &parts);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(ierr != 0, "Error in ZoltanSearchTree::ZoltanSearchTree(), call to Zoltan_LB_Partition() returned a nonzero error code.");

  for (int I = 0; I < numPoints; I++)
    partToCoordIdx[parts[I]] = I;

  Zoltan_LB_Free_Part(&impgid, &implid, &impproc, &imppart);
  Zoltan_LB_Free_Part(&gid, &lid, &procs, &parts);
}


PeridigmNS::ZoltanSearchTree::~ZoltanSearchTree() {
  Zoltan_Destroy(&zoltan);
  delete [] partToCoordIdx;
  delete [] searchParts;
}

/*
 * callback for num_objects
 */
int PeridigmNS::
ZoltanSearchTree::get_num_points(void *data,int *ierr)
{
  *ierr=ZOLTAN_OK;
  callback_data *gridData = (callback_data *)data;
  return gridData->num_points;
}

/*
 * callback for dimension
 */
int PeridigmNS::
ZoltanSearchTree::get_dimension(void *unused, int *ierr)
{
  *ierr=ZOLTAN_OK;
  return 3;
}

/*
 * callback for point ids
 */
void PeridigmNS::ZoltanSearchTree::get_point_ids(
    void *data,
    int numGids,
    int numLids,
    ZOLTAN_ID_PTR zoltanGlobalIds,
    ZOLTAN_ID_PTR zoltanLocalIds,
    int numWeights,
    float *objectWts,
    int *ierr
){
  *ierr=ZOLTAN_OK;
  callback_data *gridData = (callback_data *)data;

  /*
   * get local ids
   */
  for(int i=0;i<gridData->num_points;i++){
    zoltanGlobalIds[i]=i;
  }

}

/*
 * get point coordinates
 */
void PeridigmNS::ZoltanSearchTree::get_point_coordinates(
    void *data,
    int numGids,
    int numLids,
    int numPoints,
    ZOLTAN_ID_PTR zoltanGlobalIds,
    ZOLTAN_ID_PTR zoltanLocalIds,
    int dimension,
    double *zoltan_gridData,
    int *ierr
    ){
  /*
   * point coordinates
   */
  *ierr=ZOLTAN_OK;
  callback_data *gridData = (callback_data *)data;

  double *y=gridData->x;
  for(int i=0;i<gridData->num_points*3;i++){
    zoltan_gridData[i]=y[i];
  }

}



void PeridigmNS::ZoltanSearchTree::FindPointsWithinRadius(const double* point, double searchRadius, std::vector<int>& neighborList)
{

  /*
   * local search and therefore number of search procs = 1
   */
  int searchProcs,numSearchProcs,numSearchParts;
  Zoltan_LB_Box_PP_Assign (
      zoltan,
      *point-searchRadius,
      *(point+1)-searchRadius,
      *(point+2)-searchRadius,
      *point+searchRadius,
      *(point+1)+searchRadius,
      *(point+2)+searchRadius,
      &searchProcs,
      &numSearchProcs,
      searchParts,
      &numSearchParts);
  /*
   * filter box search to radius
   */
  double R2=searchRadius*searchRadius;
  double *X=callbackdata.x;
  for(int p=0;p<numSearchParts;p++){
    int idx=partToCoordIdx[searchParts[p]];
    double xi=X[idx*3];double yi=X[idx*3+1];double zi=X[idx*3+2];
    double dx=(xi-point[0]);
    double dy=(yi-point[1]);
    double dz=(zi-point[2]);
    if(dx*dx+dy*dy+dz*dz <=R2)
      neighborList.push_back(idx);
  }
}
