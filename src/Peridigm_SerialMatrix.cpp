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

#include <Teuchos_Exceptions.hpp>
#include <Epetra_Import.h>
#include <Epetra_SerialDenseMatrix.h>
#include "Peridigm_SerialMatrix.hpp"
#include <vector>

using namespace std;

PeridigmNS::SerialMatrix::SerialMatrix(Teuchos::RCP<Epetra_FECrsMatrix> epetraFECrsMatrix,
                                       Teuchos::RCP<const Epetra_BlockMap> epetraOverlapMap)
  : FECrsMatrix(epetraFECrsMatrix), overlapMap(epetraOverlapMap)
{
}

void PeridigmNS::SerialMatrix::addValue(int row, int col, double value)
{
  // addValue sums into the underlying Epetra_FECrsMatrix one value at a time.
  // Useful for testing, but shockingly inefficient, addValues() is prefered.

  int numRows = 1;
  int globalRowID = 3*overlapMap->GID(row/3) + row%3;
  int numCols = 1;
  int globalColID = 3*overlapMap->GID(col/3) + col%3;
  double** data = new double*[1];
  data[0] = new double[1];
  data[0][0] = value;

  int err = FECrsMatrix->SumIntoGlobalValues(numRows, &globalRowID, numCols, &globalColID, data);
  TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::SerialMatrix::addValue(), SumIntoGlobalValues() returned nonzero error code.\n");

  delete[] data[0];
  delete[] data;
}

void PeridigmNS::SerialMatrix::addValues(int numIndices, const int* indices, const double *const * values)
{
  vector<int> globalIndices(numIndices);
  vector<int> localRowIndices(numIndices);
  vector<int> localColIndices(numIndices);
  for(int i=0 ; i<numIndices ; ++i){
    int globalIndex = 3*overlapMap->GID(indices[i]/3) + indices[i]%3;
    globalIndices[i] = globalIndex;
    localRowIndices[i] = FECrsMatrix->LRID(globalIndex);
    int localColIndex = FECrsMatrix->LCID(globalIndex);
    TEST_FOR_EXCEPT_MSG(localColIndex == -1, "Error in PeridigmNS::SerialMatrix::addValues(), bad column index.");
    localColIndices[i] = localColIndex;
  }

  for(int iRow=0 ; iRow<numIndices ; ++iRow){

    // If the row is locally owned, then sum into the global tangent with Epetra_CrsMatrix::SumIntoMyValues().
    if(localRowIndices[iRow] != -1){
      int err = FECrsMatrix->SumIntoMyValues(localRowIndices[iRow], numIndices, values[iRow], &localColIndices[0]);
      TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::SerialMatrix::addValues(), SumIntoMyValues() returned nonzero error code.\n");
    }
    // If the row is not locally owned, then sum into the global tangent with Epetra_FECrsMatrix::SumIntoGlobalValues().
    // This is expensive.
    else{
      int err = FECrsMatrix->SumIntoGlobalValues(globalIndices[iRow], numIndices, values[iRow], &globalIndices[0]);
      TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::SerialMatrix::addValues(), SumIntoGlobalValues() returned nonzero error code.\n");
    }
  }
}

void PeridigmNS::SerialMatrix::putScalar(double value)
{
  FECrsMatrix->PutScalar(value);
}
