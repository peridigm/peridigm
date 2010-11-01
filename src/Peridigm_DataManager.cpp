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

#include "Peridigm_DataManager.hpp"

void PeridigmNS::DataManager::allocateData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > specs)
{
  fieldSpecs = specs;
  
  // loop over the specs and determine:
  // 1) the number of scalar, vector2d, and vector3d fields
  // 2) the FieldType for each of the data
  // 3) whether the data has one or two states
  numStatelessScalar = 0;
  numStatelessVector2D = 0;
  numStatelessVector3D = 0;
  numStatefulScalar = 0;
  numStatefulVector2D = 0;
  numStatefulVector3D = 0;
  for(unsigned int i=0; i<fieldSpecs->size() ; ++i){
    Field_NS::FieldSpec& spec = (*fieldSpecs)[i];
    if(spec.getLength() == Field_NS::FieldSpec::SCALAR){
      if(spec.getStateArchitecture() == Field_NS::FieldSpec::STATELESS)
        numStatelessScalar++;
      else
        numStatefulScalar++;
    }
    if(spec.getLength() == Field_NS::FieldSpec::VECTOR2D){
      if(spec.getStateArchitecture() == Field_NS::FieldSpec::STATELESS)
        numStatelessVector2D++;
      else
        numStatefulVector2D++;
    }
    if(spec.getLength() == Field_NS::FieldSpec::VECTOR3D){
      if(spec.getStateArchitecture() == Field_NS::FieldSpec::STATELESS)
        numStatelessVector3D++;
      else
        numStatefulVector3D++;
    }
  }

  cout << "\nDEBUGGING:" << endl;
  cout << "  numStatelessScalar " << numStatelessScalar << endl;
  cout << "  numStatefulScalar " << numStatefulScalar << endl;
  cout << "  numStatelessVector2D  " << numStatelessVector2D << endl;
  cout << "  numStatefulVector2D  " << numStatefulVector2D << endl;
  cout << "  numStatelessVector3D  " << numStatelessVector3D << endl;
  cout << "  numStatefulVector3D  " << numStatefulVector3D << endl;
  cout << endl;
  
}
