/*! \file Peridigm_State.hpp */

// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? 
// David J. Littlewood   djlittl@sandia.gov 
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

#ifndef PERIDIGM_STATE_HPP
#define PERIDIGM_STATE_HPP

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

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
  void allocateScalarData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map)
  {
    scalarData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldSpecs->size()));
    for(unsigned int i=0 ; i<fieldSpecs->size() ; ++i)
      fieldSpecToDataMap[(*fieldSpecs)[i]] = Teuchos::rcp((*scalarData)(i), false);
  }

  //! Instantiates the Epetra_MultiVector that holds vector data and creates the FieldSpec-to-Epetra_Vector mappings.
  void allocateVectorData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map)
  {
    vectorData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldSpecs->size()));
    for(unsigned int i=0 ; i<fieldSpecs->size() ; ++i)
      fieldSpecToDataMap[(*fieldSpecs)[i]] = Teuchos::rcp((*vectorData)(i), false);
  }

  //! Instantiates the Epetra_MultiVector that holds bond data and creates the FieldSpec-to-Epetra_Vector mappings.
  void allocateBondData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map)
  {
    bondData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldSpecs->size()));
    for(unsigned int i=0 ; i<fieldSpecs->size() ; ++i)
      fieldSpecToDataMap[(*fieldSpecs)[i]] = Teuchos::rcp((*bondData)(i), false);
  }

  //! Provides access to an Epetra_Vector corresponding to the given FieldSpec.
  Teuchos::RCP<Epetra_Vector> getData(Field_NS::FieldSpec fieldSpec){

    // search for the data
    std::map< Field_NS::FieldSpec, Teuchos::RCP<Epetra_Vector> >::iterator lb = fieldSpecToDataMap.lower_bound(fieldSpec);
    // if the key does not exist, throw an exception
    bool keyExists = ( lb != fieldSpecToDataMap.end() && !(fieldSpecToDataMap.key_comp()(fieldSpec, lb->first)) );
    TEST_FOR_EXCEPTION(!keyExists, Teuchos::RangeError, 
                       "Error in PeridigmNS::State::getData(), key does not exist!");
    return lb->second;
  }

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
