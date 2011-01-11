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
  DataManager(const DataManager& dataManager){
    // \todo Write me.
  }
  ~DataManager(){}

  void setMaps(Teuchos::RCP<const Epetra_BlockMap> scalarMap_,
               Teuchos::RCP<const Epetra_BlockMap> vector3DMap_,
               Teuchos::RCP<const Epetra_BlockMap> bondMap_){
    scalarMap = scalarMap_;
    vector3DMap = vector3DMap_;
    bondMap = bondMap_;
  }

  void allocateData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > specs);

  void rebalance(Teuchos::RCP<const Epetra_BlockMap> rebalancedScalarMap,
                 Teuchos::RCP<const Epetra_BlockMap> rebalancedVector3DMap,
                 Teuchos::RCP<const Epetra_BlockMap> rebalancedBondMap);

  Teuchos::RCP<Epetra_Vector> getData(Field_NS::FieldSpec fieldSpec,
                                      Field_NS::FieldSpec::FieldStep fieldStep);

  void updateState(){
    stateN.swap(stateNP1);
  }

protected:

  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessScalarFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessVector3DFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessBondFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulScalarFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulVector3DFieldSpecs;
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulBondFieldSpecs;
  Teuchos::RCP<const Epetra_BlockMap> scalarMap;
  Teuchos::RCP<const Epetra_BlockMap> vector3DMap;
  Teuchos::RCP<const Epetra_BlockMap> bondMap;
  Teuchos::RCP<State> stateN;
  Teuchos::RCP<State> stateNP1;
  Teuchos::RCP<State> stateNONE;
};

}

#endif // PERIDIGM_DATAMANAGER_HPP
