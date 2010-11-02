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

class State {

public:

  State(){}
  State(const State& state){}
  ~State(){}

  void allocateScalarData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map)
  {
    scalarData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldSpecs->size()));
    for(unsigned int i=0 ; i<fieldSpecs->size() ; ++i)
      fieldSpecToDataMap[(*fieldSpecs)[i]] = Teuchos::rcp((*scalarData)(i), false);
  }

  void allocateVector2DData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map)
  {
    vector2DData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldSpecs->size()));
    for(unsigned int i=0 ; i<fieldSpecs->size() ; ++i)
      fieldSpecToDataMap[(*fieldSpecs)[i]] = Teuchos::rcp((*vector2DData)(i), false);
  }

  void allocateVector3DData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map)
  {
    vector3DData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldSpecs->size()));
    for(unsigned int i=0 ; i<fieldSpecs->size() ; ++i)
      fieldSpecToDataMap[(*fieldSpecs)[i]] = Teuchos::rcp((*vector3DData)(i), false);
  }

  void allocateBondData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs, Teuchos::RCP<const Epetra_BlockMap> map)
  {
    bondData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldSpecs->size()));
    for(unsigned int i=0 ; i<fieldSpecs->size() ; ++i)
      fieldSpecToDataMap[(*fieldSpecs)[i]] = Teuchos::rcp((*bondData)(i), false);
  }

  Teuchos::RCP<Epetra_Vector> getData(Field_NS::FieldSpec fieldSpec){
    // search for the data
    std::map< Field_NS::FieldSpec, Teuchos::RCP<Epetra_Vector> >::iterator lb = fieldSpecToDataMap.lower_bound(fieldSpec);
    // if the key does not exist, throw an exception
    bool keyExists = ( lb != fieldSpecToDataMap.end() && !(fieldSpecToDataMap.key_comp()(fieldSpec, lb->first)) );
    TEST_FOR_EXCEPTION(!keyExists, Teuchos::RangeError, 
                       "Error in PeridigmNS::State::getData(), key does not exist!");
    return lb->second;
  }

protected:

  Teuchos::RCP<Epetra_MultiVector> scalarData;
  Teuchos::RCP<Epetra_MultiVector> vector2DData;
  Teuchos::RCP<Epetra_MultiVector> vector3DData;
  Teuchos::RCP<Epetra_MultiVector> bondData;
  std::map< Field_NS::FieldSpec, Teuchos::RCP<Epetra_Vector> > fieldSpecToDataMap;
};

}

#endif // PERIDIGM_STATE_HPP
