/*! \file Peridigm_DataManager.cpp */

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

#include <Teuchos_Exceptions.hpp>
#include "Peridigm_DataManager.hpp"

void PeridigmNS::DataManager::allocateData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > specs)
{
  fieldSpecs = specs;
  
  // loop over the specs and determine:
  // 1) the number of scalar, vector2d, and vector3d fields
  // 2) the FieldType for each of the data
  // 3) whether the data has one or two states
  statelessScalarFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statelessVector2DFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statelessVector3DFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statefulScalarFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statefulVector2DFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statefulVector3DFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  for(unsigned int i=0; i<fieldSpecs->size() ; ++i){
    Field_NS::FieldSpec& spec = (*fieldSpecs)[i];
    if(spec.getLength() == Field_NS::FieldSpec::SCALAR){
      if(spec.getStateArchitecture() == Field_NS::FieldSpec::STATELESS)
        statelessScalarFieldSpecs->push_back(spec);
      else
        statefulScalarFieldSpecs->push_back(spec);
    }
    if(spec.getLength() == Field_NS::FieldSpec::VECTOR2D){
      if(spec.getStateArchitecture() == Field_NS::FieldSpec::STATELESS)
        statelessVector2DFieldSpecs->push_back(spec);
      else
        statefulVector2DFieldSpecs->push_back(spec);
    }
    if(spec.getLength() == Field_NS::FieldSpec::VECTOR3D){
      if(spec.getStateArchitecture() == Field_NS::FieldSpec::STATELESS)
        statelessVector3DFieldSpecs->push_back(spec);
      else
        statefulVector3DFieldSpecs->push_back(spec);
    }
  }

  cout << "\nDEBUGGING:" << endl;
  cout << "  numStatelessScalar " << statelessScalarFieldSpecs->size() << endl;
  cout << "  numStatefulScalar " << statefulScalarFieldSpecs->size() << endl;
  cout << "  numStatelessVector2D  " << statelessVector2DFieldSpecs->size() << endl;
  cout << "  numStatefulVector2D  " << statefulVector2DFieldSpecs->size() << endl;
  cout << "  numStatelessVector3D  " << statelessVector3DFieldSpecs->size() << endl;
  cout << "  numStatefulVector3D  " << statefulVector3DFieldSpecs->size() << endl;
  cout << endl;

  // make sure maps exist before trying to create states
  if(statelessScalarFieldSpecs->size() + statefulScalarFieldSpecs->size() > 0)
    TEST_FOR_EXCEPTION(scalarMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate scalar data with no scalar data map (forget setScalarMap()?).");
  if(statelessVector2DFieldSpecs->size() + statefulVector2DFieldSpecs->size() > 0)
    TEST_FOR_EXCEPTION(vector2DMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate vector2D data with no scalar data map (forget setVector2DMap()?).");
  if(statelessVector3DFieldSpecs->size() != 0 + statefulVector3DFieldSpecs->size() > 0)
    TEST_FOR_EXCEPTION(vector3DMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate vector3D data with no scalar data map (forget setVector3DMap()?).");

  // create the states
  if(statelessScalarFieldSpecs->size() + statelessVector2DFieldSpecs->size() + statelessVector3DFieldSpecs->size() > 0){
    stateNONE = Teuchos::rcp(new State);
    if(statelessScalarFieldSpecs->size() > 0)
      stateNONE->allocateScalarData(statelessScalarFieldSpecs, scalarMap);
    if(statelessVector2DFieldSpecs->size() > 0)
      stateNONE->allocateVector2DData(statelessVector2DFieldSpecs, vector2DMap);
    if(statelessVector3DFieldSpecs->size() > 0)
      stateNONE->allocateVector3DData(statelessVector3DFieldSpecs, vector3DMap);
  }
  if(statefulScalarFieldSpecs->size() + statefulVector2DFieldSpecs->size() + statefulVector3DFieldSpecs->size() > 0){
    stateN = Teuchos::rcp(new State);
    stateNP1 = Teuchos::rcp(new State);
    if(statefulScalarFieldSpecs->size() > 0){
      stateN->allocateScalarData(statefulScalarFieldSpecs, scalarMap);
      stateNP1->allocateScalarData(statefulScalarFieldSpecs, scalarMap);
    }
    if(statefulVector2DFieldSpecs->size() > 0){
      stateN->allocateVector2DData(statefulVector2DFieldSpecs, vector2DMap);
      stateNP1->allocateVector2DData(statefulVector2DFieldSpecs, vector2DMap);
    }
    if(statefulVector3DFieldSpecs->size() > 0){
      stateN->allocateVector3DData(statefulVector3DFieldSpecs, vector3DMap);
      stateNP1->allocateVector3DData(statefulVector3DFieldSpecs, vector3DMap);
    }
  }
}

Teuchos::RCP<Epetra_Vector> PeridigmNS::DataManager::getData(Field_NS::FieldSpec fieldSpec, Field_NS::FieldSpec::FieldStep fieldStep)
{
  Teuchos::RCP<Epetra_Vector> data;
  if(fieldStep == Field_NS::FieldSpec::STEP_NONE){
    data = stateNONE->getData(fieldSpec);
  }
  else if(fieldStep == Field_NS::FieldSpec::STEP_N){
    data = stateN->getData(fieldSpec);
  }
  else if(fieldStep == Field_NS::FieldSpec::STEP_NP1){
    data = stateNP1->getData(fieldSpec);
  }
  else{
    TEST_FOR_EXCEPTION(false, Teuchos::RangeError, 
                       "PeridigmNS::DataManager::getData, invalid FieldStep!");
  }
  return data;
}
