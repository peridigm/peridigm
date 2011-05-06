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
#include "Field.h"

namespace PeridigmNS {

/*! \brief Container class for scalar, vector, and bond data.
 *
 * The State class is a container class for scalar, vector, and bond data.  Data is stored in Epetra_MultiVector objects.
 * Individual Epetra_Vectors corresponding to a specific FieldSpec are accessed through the getData function.  The
 * State class was designed for use within a DataManger.
 */
class State {

public:

  //! Constructor.
  State(){}

  //! Copy constructor.
  State(const State& state){}

  //! Destructor.
  ~State(){}

  //! Instantiates the Epetra_MultiVector that holds scalar data and creates the FieldSpec-to-Epetra_Vector mappings.
  void allocateScalarData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map);

  //! Instantiates the Epetra_MultiVector that holds vector data and creates the FieldSpec-to-Epetra_Vector mappings.
  void allocateVectorData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map);

  //! Instantiates the Epetra_MultiVector that holds bond data and creates the FieldSpec-to-Epetra_Vector mappings.
  void allocateBondData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map);

  //! Provides access to an Epetra_Vector corresponding to the given FieldSpec.
  Teuchos::RCP<Epetra_Vector> getData(Field_NS::FieldSpec fieldSpec);

  //! @name Epetra_MultiVector accessor functions
  //@{
  //! Accessor for the scalar Epetra_MultiVector.
  Teuchos::RCP<Epetra_MultiVector> getScalarMultiVector() { return scalarData; }
  //! Accessor for the vector Epetra_MultiVector.
  Teuchos::RCP<Epetra_MultiVector> getVectorMultiVector() { return vectorData; }
  //! Accessor for the bond-data Epetra_MultiVector.
  Teuchos::RCP<Epetra_MultiVector> getBondMultiVector() { return bondData; }
  //@}

protected:

  //! @name Epetra_MultiVectors
  //@{
  //! Epetra_MultiVector for scalar data.
  Teuchos::RCP<Epetra_MultiVector> scalarData;
  //! Epetra_MultiVector for vector data.
  Teuchos::RCP<Epetra_MultiVector> vectorData;
  //! Epetra_MultiVector for bond data.
  Teuchos::RCP<Epetra_MultiVector> bondData;
  //@}

  //! Map that associates a FieldSpec with an individual Epetra_Vector contained within one of the Epetra_MultiVectors.
  std::map< Field_NS::FieldSpec, Teuchos::RCP<Epetra_Vector> > fieldSpecToDataMap;
};

}

#endif // PERIDIGM_STATE_HPP
