//*! \file Peridigm_ComputeManager.cpp */
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

#include <vector>
#include <string>
#include <iostream>

#include "Peridigm_ComputeManager.hpp"
#include "compute/Peridigm_Compute.hpp"
#include "compute/Peridigm_Compute_ACCEL3D.hpp"
#include "Peridigm_DataManager.hpp"
#include <Field.h>

PeridigmNS::ComputeManager::ComputeManager( Teuchos::RCP<Teuchos::ParameterList>& params) {

  Teuchos::RCP<Compute> compute;

  // No input to validate; no computes requested
  if (params == Teuchos::null) return;

  if (params->isParameter("Acceleration")) {
    compute = Teuchos::rcp(new PeridigmNS::Compute_ACCEL3D );
    computeObjects.push_back( Teuchos::rcp_implicit_cast<Compute>(compute) );
  }
}

Teuchos::ParameterList PeridigmNS::ComputeManager::getValidParameterList() {
  Teuchos::ParameterList validParameterList("Output");
  return validParameterList;
}

std::vector<Field_NS::FieldSpec> PeridigmNS::ComputeManager::getFieldSpecs() {

  std::vector<Field_NS::FieldSpec> myFieldSpecs;

  // Loop over all compute objects, collect the field specs they compute
  for (unsigned int i=0; i < computeObjects.size(); i++) {
    Teuchos::RCP<const PeridigmNS::Compute> compute = computeObjects[i];
    std::vector<Field_NS::FieldSpec> computeFieldSpecs = compute->getFieldSpecs();
    myFieldSpecs.insert(myFieldSpecs.end(), computeFieldSpecs.begin(), computeFieldSpecs.end());
  }

  // remove duplicates
  std::sort(myFieldSpecs.begin(), myFieldSpecs.end());
  std::unique(myFieldSpecs.begin(), myFieldSpecs.end());

  return myFieldSpecs;
}

PeridigmNS::ComputeManager::~ComputeManager() {
}

void PeridigmNS::ComputeManager::compute(Teuchos::RCP<PeridigmNS::DataManager>& dataManager) {
}
