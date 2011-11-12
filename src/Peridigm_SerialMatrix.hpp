/*! \file Peridigm_SerialMatrix.hpp */

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

#ifndef PERIDIGM_DENSEMATRIX_HPP
#define PERIDIGM_DENSEMATRIX_HPP

#include <vector>
#include <map>
#include <iostream>
#include <Teuchos_RCP.hpp>
#include <Epetra_FECrsMatrix.h>

namespace PeridigmNS {

/*! \brief A dense serial matrix.
 */
class SerialMatrix {

public:

  SerialMatrix(Teuchos::RCP<Epetra_FECrsMatrix> epetraFECrsMatrix,
               Teuchos::RCP<const Epetra_BlockMap> epetraOverlapMap);

  //! Destructor.
  ~SerialMatrix(){}

  //! Add data at given location, indexed by local ID
  void addValue(int row, int col, double value);

  //! Add block of data at given locations, indexed by local ID
  void addValues(int numIndicies, const int* indicies, const double *const * values);

  //! Set all entries to given scalar
  void putScalar(double value);

  //! Return ref-count pointer to the FECrsMatrix
  Teuchos::RCP<const Epetra_FECrsMatrix> getFECrsMatrix() { return FECrsMatrix; }

protected:

  Teuchos::RCP<Epetra_FECrsMatrix> FECrsMatrix;
  Teuchos::RCP<const Epetra_BlockMap> overlapMap;

private:

  //! Private to prohibit use.
  SerialMatrix() {}
  SerialMatrix(const SerialMatrix& serialMatrix){}
};

}

#endif // PERIDIGM_DENSEMATRIX_HPP
