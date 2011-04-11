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
#include "Peridigm_SerialMatrix.hpp"
#include <sstream>

PeridigmNS::SerialMatrix::SerialMatrix(Teuchos::RCP<Epetra_FECrsMatrix> epetraFECrsMatrix,
                                       Teuchos::RCP<const Epetra_BlockMap> epetraOverlapMap)
  : FECrsMatrix(epetraFECrsMatrix), overlapMap(epetraOverlapMap)
{
}

void PeridigmNS::SerialMatrix::addValue(int row, int col, double value)
{
  // what we really want to do here is cache a block of values and their
  // row and column ids.  Then when the block is full, sum into the crs matrix.

  // the initial implementation is just an interface to the underlying Epetra_FECrsMatrix,
  // filling one value at a time (very inefficient).

  int numRows = 1;
  int globalRowID = 3*overlapMap->GID(row/3) + row%3;
  int numCols = 1;
  int globalColID = 3*overlapMap->GID(col/3) + col%3;
  double** data = new double*[1];
  data[0] = new double[1];
  data[0][0] = value;

  FECrsMatrix->SumIntoGlobalValues(numRows, &globalRowID, numCols, &globalColID, data);

  delete[] data[0];
  delete[] data;
}

void PeridigmNS::SerialMatrix::putScalar(double value)
{
  FECrsMatrix->PutScalar(value);
}
