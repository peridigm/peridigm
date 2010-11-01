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

  DataManager() : numStatelessScalar(0), numStatelessVector2D(0), numStatelessVector3D(0),
                  numStatefulScalar(0), numStatefulVector2D(0), numStatefulVector3D(0) {}
  DataManager(const DataManager& dataManager){}
  ~DataManager(){}

  void setScalarMap(Teuchos::RCP<Epetra_BlockMap> map) { scalarMap = map; }
  void setVector2DMap(Teuchos::RCP<Epetra_BlockMap> map) { vector2DMap = map; }
  void setVector3DMap(Teuchos::RCP<Epetra_BlockMap> map) { vector3DMap = map; }
  void allocateData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > specs);

protected:

  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs;
  int numStatelessScalar, numStatelessVector2D, numStatelessVector3D;
  int numStatefulScalar, numStatefulVector2D, numStatefulVector3D;
  Teuchos::RCP<Epetra_BlockMap> scalarMap;
  Teuchos::RCP<Epetra_BlockMap> vector2DMap;
  Teuchos::RCP<Epetra_BlockMap> vector3DMap;
  Teuchos::RCP<State> stateN;
  Teuchos::RCP<State> stateNP1;
  Teuchos::RCP<State> stateNONE;
  std::map< Field_NS::FieldSpec, Teuchos::RCP<Epetra_Vector> > fieldSpecToDataMap;
};

}

#endif // PERIDIGM_DATAMANAGER_HPP
