/*! \file Peridigm_State.hpp */

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

#ifndef PERIDIGM_STATE_HPP
#define PERIDIGM_STATE_HPP

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include "Peridigm_Field.hpp"
#include <vector>

namespace PeridigmNS {

/*! \brief Container class for scalar point data, vector point data, and scalar bond data.
 *
 * The State class is a container class for storing fields of various length (scalar, vector, bond, etc.).
 * Data is stored in Epetra_MultiVector objects.  Individual Epetra_Vectors corresponding to a specific
 * field id are accessed through the getData function.  The State class was designed for use within a
 * DataManger.
 */
class State {

public:

  //! Constructor.
  State() : 
    maxPointDataElementSize(9),
    numFieldIds(0),
    pointData(std::vector< Teuchos::RCP<Epetra_MultiVector> >(maxPointDataElementSize)) {}

  //! Copy constructor.
  State(const State& state) {}
 
  //! Destructor.
  ~State() {}

  //! @name Memory allocation functions
  //@{

  /** \brief Allocates underlying Epetra_Multivectors for point-wise data.
  **
  **  The length argument denotes data size (e.g., 1 = scalar, 3 = vector, etc.).  This function may be called only
  **  once for a given data size.  The fieldIds and map provided must be consistent with the given length.
  **/
  void allocatePointData(PeridigmField::Length length, std::vector<int> fieldIds, Teuchos::RCP<const Epetra_BlockMap> map);

  //! Allocates underlying Epetra_Multivector for bond data; only scalar bond data is supported.
  void allocateBondData(std::vector<int> fieldIds, Teuchos::RCP<const Epetra_BlockMap> map);

  //@}

  //! Return the maximum allowable element size for point data.
  int getMaxPointDataElementSize(){ return maxPointDataElementSize; }

  //! @name Accessor functions for underlying Epetra_Multivectors.
  //@{

  //! Accessor for underlying point-wise data Epetra_MultiVectors.
  Teuchos::RCP<Epetra_MultiVector> getPointMultiVector(PeridigmField::Length length) { return pointData[PeridigmField::variableDimension(length)-1]; }

  //! Accessor for underlying point-wise data Epetra_MultiVectors.
  Teuchos::RCP<Epetra_MultiVector> getPointMultiVector(int index) { return pointData[index]; }

  //! Accessor for underlying bond data Epetra_MultiVectors.
  Teuchos::RCP<Epetra_MultiVector> getBondMultiVector() { return bondData; }

  //@}

  //! Returns the list of field ids of the given relation and length; if no arguments are given, returns complete list of field ids.
  std::vector<int> getFieldIds(PeridigmField::Relation relation = PeridigmField::UNDEFINED_RELATION,
                               PeridigmField::Length length = PeridigmField::UNDEFINED_LENGTH);

  //! Query the existence of a field id.
  bool hasData(int fieldId);

  //! Provides access to an Epetra_Vector corresponding to the given field id.
  Teuchos::RCP<Epetra_Vector> getData(int fieldId);

  //! Copies data from a different state object based on global IDs; functions only if all the local IDs in the target map exist in and are locally owned in the source map.
  void copyLocallyOwnedDataFromState(Teuchos::RCP<PeridigmNS::State> source);

private:

  //! Maximum set of an element in the pointData Epetra_MultiVector.
  int maxPointDataElementSize;

protected:

  //! Maximum value of field ids.
  int numFieldIds;

  //! Epetra_MultiVectors for point data.
  std::vector< Teuchos::RCP<Epetra_MultiVector> > pointData;

  //! Epetra_MultiVector for bond data.
  Teuchos::RCP<Epetra_MultiVector> bondData;

  //! Map that associates a field id with an individual Epetra_Vector contained within one of the Epetra_MultiVectors.
  std::map< int, Teuchos::RCP<Epetra_Vector> > fieldIdToDataMap;

  //! Vector that associates a field id with an individual Epetra_Vector contained within one of the Epetra_MultiVectors.
  //  The vector is used in addition to the map because access is faster in getData() and this was shown to be a bottleneck.
  std::vector< Teuchos::RCP<Epetra_Vector> > fieldIdToDataVector;

  //! Utility for copying locally-owned data from one multivector to another; succeeds only if all the local IDs in the target map exist in and are locally owned in the source map.
  void copyLocallyOwnedMultiVectorData(Epetra_MultiVector& source, Epetra_MultiVector& target);
};

}

#endif // PERIDIGM_STATE_HPP
