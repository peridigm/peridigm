/*! \file Peridigm_DataManager.hpp */

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

#ifndef PERIDIGM_DATAMANAGER_HPP
#define PERIDIGM_DATAMANAGER_HPP

#include "Field.h"
#include "Peridigm_State.hpp"

namespace PeridigmNS {

class DataManager {

public:

  DataManager(){}
  DataManager(const DataManager& dataManager){}
  ~DataManager(){}

  void setScalarMap(Teuchos::RCP<const Epetra_BlockMap> map) { scalarMap = map; }
  void setVector2DMap(Teuchos::RCP<const Epetra_BlockMap> map) { vector2DMap = map; }
  void setVector3DMap(Teuchos::RCP<const Epetra_BlockMap> map) { vector3DMap = map; }
  void setBondMap(Teuchos::RCP<const Epetra_BlockMap> map) { bondMap = map; }
  void allocateData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > specs);
  Teuchos::RCP<Epetra_Vector> getData(Field_NS::FieldSpec fieldSpec, Field_NS::FieldSpec::FieldStep fieldStep);

  void updateState(){
    Teuchos::RCP<State> temp = stateN;
    stateN = stateNP1;
    stateNP1 = temp;
  }

protected:

  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessScalarFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessVector2DFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessVector3DFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessBondFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulScalarFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulVector2DFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulVector3DFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulBondFieldSpecs;
  Teuchos::RCP<const Epetra_BlockMap> scalarMap;
  Teuchos::RCP<const Epetra_BlockMap> vector2DMap;
  Teuchos::RCP<const Epetra_BlockMap> vector3DMap;
  Teuchos::RCP<const Epetra_BlockMap> bondMap;
  Teuchos::RCP<State> stateN;
  Teuchos::RCP<State> stateNP1;
  Teuchos::RCP<State> stateNONE;
};

}

#endif // PERIDIGM_DATAMANAGER_HPP
