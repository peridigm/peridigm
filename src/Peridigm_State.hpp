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
#include "mesh_output/vtk/Field.h"

namespace PeridigmNS {

/*! \brief Container class for scalar point data, vector point data, and scalar bond data.
 *
 * The State class is a container class for scalar point data, vector point data, and scalar bond data.
 * Data is stored in Epetra_MultiVector objects.  Individual Epetra_Vectors corresponding to a specific
 * FieldSpec are accessed through the getData function.  The State class was designed for use within a
 * DataManger.
 */
class State {

public:

  //! Constructor.
  State(){}

  //! Copy constructor.
  State(const State& state){}

  //! Destructor.
  ~State(){}

  //! @name Memory allocation functions
  //@{

  //! Instantiates the Epetra_MultiVector that holds scalar point data and creates the FieldSpec-to-Epetra_Vector mappings.
  void allocateScalarPointData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map);

  //! Instantiates the Epetra_MultiVector that holds vector point data and creates the FieldSpec-to-Epetra_Vector mappings.
  void allocateVectorPointData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map);

  //! Instantiates the Epetra_MultiVector that holds scalar bond data and creates the FieldSpec-to-Epetra_Vector mappings.
  void allocateScalarBondData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map);

  //@}

  //! @name Epetra_MultiVector accessor functions
  //@{

  //! Accessor for the scalar point Epetra_MultiVector.
  Teuchos::RCP<Epetra_MultiVector> getScalarPointMultiVector() { return scalarPointData; }

  //! Accessor for the vector point Epetra_MultiVector.
  Teuchos::RCP<Epetra_MultiVector> getVectorPointMultiVector() { return vectorPointData; }

  //! Accessor for the scalar bond Epetra_MultiVector.
  Teuchos::RCP<Epetra_MultiVector> getScalarBondMultiVector() { return scalarBondData; }

  //@}

  //! Returns the list of field specs of the given relation (POINT, BOND) and length (SCALAR, VECTOR3D); if no argument is given, returns complete list of field specs.
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > getFieldSpecs(Teuchos::RCP<Field_ENUM::Relation> relation = Teuchos::RCP<Field_ENUM::Relation>(),
		                                                         Teuchos::RCP<Field_ENUM::Length> length = Teuchos::RCP<Field_ENUM::Length>());

  //! Query the existence of a fieldSpec.
  bool hasData(Field_NS::FieldSpec fieldSpec);

  //! Provides access to an Epetra_Vector corresponding to the given FieldSpec.
  Teuchos::RCP<Epetra_Vector> getData(Field_NS::FieldSpec fieldSpec);

  //! Copies data from a different state object based on global IDs; functions only if all the local IDs in the target map exist in and are locally owned in the source map.
  void copyLocallyOwnedDataFromState(Teuchos::RCP<PeridigmNS::State> source);

protected:

  //! @name Epetra_MultiVectors
  //@{
  //! Epetra_MultiVector for scalar point data.
  Teuchos::RCP<Epetra_MultiVector> scalarPointData;
  //! Epetra_MultiVector for vector point data.
  Teuchos::RCP<Epetra_MultiVector> vectorPointData;
  //! Epetra_MultiVector for scalar bond data.
  Teuchos::RCP<Epetra_MultiVector> scalarBondData;
  //@}

  //! Map that associates a FieldSpec with an individual Epetra_Vector contained within one of the Epetra_MultiVectors.
  std::map< Field_NS::FieldSpec, Teuchos::RCP<Epetra_Vector> > fieldSpecToDataMap;

  //! Utility for copying locally-owned data from one multivector to another; succeeds only if all the local IDs in the target map exist in and are locally owned in the source map.
  void copyLocallyOwnedMultiVectorData(Epetra_MultiVector& source, Epetra_MultiVector& target);
};

}

#endif // PERIDIGM_STATE_HPP
