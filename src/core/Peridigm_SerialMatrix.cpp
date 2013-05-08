/*! \file Peridigm_SerialMatrix.cpp */

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

#include <vector>

#include <Epetra_Import.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_Exceptions.hpp>
#include <Teuchos_SerialDenseVector.hpp>

#include "Peridigm_SerialMatrix.hpp"

using namespace std;

PeridigmNS::SerialMatrix::SerialMatrix(Teuchos::RCP<Epetra_FECrsMatrix> epetraFECrsMatrix)
  : FECrsMatrix(epetraFECrsMatrix)
{
}

void PeridigmNS::SerialMatrix::addValue(int globalRow, int globalCol, double value)
{
  // addValue sums into the underlying Epetra_FECrsMatrix one value at a time.
  // Useful for testing, but shockingly inefficient, addValues() is prefered.

  int numRows = 1;
  int numCols = 1;
  double** data = new double*[1];
  data[0] = new double[1];
  data[0][0] = value;

  int err = FECrsMatrix->SumIntoGlobalValues(numRows, &globalRow, numCols, &globalCol, data);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::SerialMatrix::addValue(), SumIntoGlobalValues() returned nonzero error code.\n");

  delete[] data[0];
  delete[] data;
}

void PeridigmNS::SerialMatrix::addValues(int numIndices, const int* globalIndices, const double *const * values)
{
  vector<int> localRowIndices(numIndices);
  vector<int> localColIndices(numIndices);
  for(int i=0 ; i<numIndices ; ++i){
    localRowIndices[i] = FECrsMatrix->LRID(globalIndices[i]);
    int localColIndex = FECrsMatrix->LCID(globalIndices[i]);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(localColIndex == -1, "Error in PeridigmNS::SerialMatrix::addValues(), bad column index.");
    localColIndices[i] = localColIndex;
  }

  for(int iRow=0 ; iRow<numIndices ; ++iRow){

    // If the row is locally owned, then sum into the global tangent with Epetra_CrsMatrix::SumIntoMyValues().
    if(localRowIndices[iRow] != -1){
      int err = FECrsMatrix->SumIntoMyValues(localRowIndices[iRow], numIndices, values[iRow], &localColIndices[0]);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::SerialMatrix::addValues(), SumIntoMyValues() returned nonzero error code.\n");
    }
    // If the row is not locally owned, then sum into the global tangent with Epetra_FECrsMatrix::SumIntoGlobalValues().
    // This is expensive.
    else{
      int err = FECrsMatrix->SumIntoGlobalValues(globalIndices[iRow], numIndices, values[iRow], &globalIndices[0]);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::SerialMatrix::addValues(), SumIntoGlobalValues() returned nonzero error code.\n");
    }
  }
}

// This is like the SerialMatrix::addValues routine above, but inserts only the block diagonal values and filters out the rest
void PeridigmNS::SerialMatrix::addBlockDiagonalValues(int numIndices, const int* globalIndices, const double *const * values)
{

  // Local row and column indices for each global index
  vector<int> localRowIndices(numIndices);
  vector<int> localColIndices(numIndices);

  // Inverse map that gives index value into globalIndices array for each global index value
  std::map<int,int> inverseMap;

  // Assign local indices and inverse map
  for(int i=0 ; i<numIndices ; ++i){
    localRowIndices[i] = FECrsMatrix->LRID(globalIndices[i]);
    int localColIndex = FECrsMatrix->LCID(globalIndices[i]);
    // Will be receiving data for columns that we will not fill, so don't check that all column data is locally owned. 
    localColIndices[i] = localColIndex;
    inverseMap[globalIndices[i]] = i;
  }

  // Scratch space for extracting the three nonzeros per row to fill
  int blockDiagonalNumIndices = 3;
  Teuchos::SerialDenseVector<int,int> blockDiagonalLocalColIndices(blockDiagonalNumIndices);
  Teuchos::SerialDenseVector<int,int> blockDiagonalGlobalIndices(blockDiagonalNumIndices);
  Teuchos::SerialDenseVector<int,double> blockDiagonalValues(blockDiagonalNumIndices);

  for(int iRow=0 ; iRow<numIndices ; ++iRow){
 
    // Determine which global element iRow belongs to
    int elem = globalIndices[iRow] / 3;
    // Determine global indices of DOFs for this element
    int e1 = 3*elem + 0;
    int e2 = 3*elem + 1;
    int e3 = 3*elem + 2;

    // Be sure entries exist in inverseMap
    TEUCHOS_TEST_FOR_EXCEPT_MSG(inverseMap.count(e1)<=0, "Error in PeridigmNS::SerialMatrix::addBlockDiagonalValues(), bad index.");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(inverseMap.count(e2)<=0, "Error in PeridigmNS::SerialMatrix::addBlockDiagonalValues(), bad index.");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(inverseMap.count(e3)<=0, "Error in PeridigmNS::SerialMatrix::addBlockDiagonalValues(), bad index.");

    // Store only the three block diagonal nonzeros to fill
    int idx1 = inverseMap[e1];
    blockDiagonalLocalColIndices[0] = localColIndices[idx1];
    blockDiagonalValues[0]          = values[iRow][idx1];
    int idx2 = inverseMap[e2];
    blockDiagonalLocalColIndices[1] = localColIndices[idx2];
    blockDiagonalValues[1]          = values[iRow][idx2];
    int idx3 = inverseMap[e3];
    blockDiagonalLocalColIndices[2] = localColIndices[idx3];
    blockDiagonalValues[2]          = values[iRow][idx3];
    // Store global indices in case row not locally owned
    blockDiagonalGlobalIndices[0] = e1;
    blockDiagonalGlobalIndices[1] = e2;
    blockDiagonalGlobalIndices[2] = e3;
     
    // If the row is locally owned, then sum into the global tangent with Epetra_CrsMatrix::SumIntoMyValues().
    if(localRowIndices[iRow] != -1){
      int err = FECrsMatrix->SumIntoMyValues(localRowIndices[iRow], blockDiagonalNumIndices, &blockDiagonalValues[0], &blockDiagonalLocalColIndices[0]);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::SerialMatrix::addBlockDiagonalValues(), SumIntoMyValues() returned nonzero error code.\n");
    }
    // If the row is not locally owned, then sum into the global tangent with Epetra_FECrsMatrix::SumIntoGlobalValues().
    // This is expensive.
    else{
      int err = FECrsMatrix->SumIntoGlobalValues(globalIndices[iRow], blockDiagonalNumIndices, &blockDiagonalValues[0], &blockDiagonalGlobalIndices[0]);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::SerialMatrix::addBlockDiagonalValues(), SumIntoGlobalValues() returned nonzero error code.\n");
    }
  }
}

void PeridigmNS::SerialMatrix::putScalar(double value)
{
  FECrsMatrix->PutScalar(value);
}
