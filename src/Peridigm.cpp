/*! \file Peridigm.cpp
 *
 * File containing main class for Peridigm: A parallel, multi-physics,
 * peridynamics simulation code.
 */

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

#include <iostream>
#include <vector>

#include <Teuchos_VerboseObject.hpp>

#include "Peridigm.hpp"
#include "Peridigm_DiscretizationFactory.hpp"
#include "contact/Peridigm_ContactModel.hpp"
#include "contact/Peridigm_ShortRangeForceContactModel.hpp"
#include "materials/Peridigm_LinearElasticIsotropicMaterial.hpp"
#include "materials/Peridigm_IsotropicElasticPlasticMaterial.hpp"

Peridigm::Peridigm::Peridigm(const Teuchos::RCP<const Epetra_Comm>& comm,
                   const Teuchos::RCP<Teuchos::ParameterList>& params)
{

  peridigmComm = comm;
  peridigmParams = params;

  bondData = NULL;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // Read mesh from disk or generate using geometric primatives.
  // All maps are generated here
  initializeDiscretization();

  // Setup contact
  initializeContact();


}

void Peridigm::Peridigm::initializeDiscretization() {

  // Extract problem parameters sublist
  Teuchos::RCP<Teuchos::ParameterList> problemParams = Teuchos::rcp(&(peridigmParams->sublist("Problem")),false);

  // Extract discretization parameters sublist
  Teuchos::RCP<Teuchos::ParameterList> discParams = Teuchos::rcp(&(problemParams->sublist("Discretization")), false);

  // Create discretization object
  // The discretization object creates:
  // 1) The epetra maps
  // 2) The vector of initial positions
  DiscretizationFactory discFactory(discParams);
  Teuchos::RCP<AbstractDiscretization> peridigmDisc = discFactory.create(peridigmComm);

  // oneDimensionalMap
  // used for cell volumes and scalar constitutive data
  oneDimensionalMap = peridigmDisc->getOneDimensionalMap(); 

  // oneDimensionalOverlapMap
  // used for cell volumes and scalar constitutive data
  // includes ghosts
  oneDimensionalOverlapMap = peridigmDisc->getOneDimensionalOverlapMap();

  // threeDimensionalMap
  // used for positions, displacements, velocities and vector constitutive data
// MLP
//  threeDimensionalMap = peridigmDisc->getThreeDimensionalMap();

  // threeDimensionalOverlapMap
  // used for positions, displacements, velocities and vector constitutive data
  // includes ghosts
// MLP
//  threeDimensionalOverlapMap = peridigmDisc->getThreeDimensionalOverlapMap();

  // bondConstitutiveDataMap
  // a non-overlapping map used for storing constitutive data on bonds
  bondMap = peridigmDisc->getBondMap();

}

void Peridigm::Peridigm::initializeContact() {

  // Extract problem parameters sublist
  Teuchos::RCP<Teuchos::ParameterList> problemParams = Teuchos::rcp(&(peridigmParams->sublist("Problem")),false);

  // Assume no contact
  computeContact = false;
  contactSearchRadius = 0.0;
  contactSearchFrequency = 0;


  // Set up contact, if requested by user
  if(problemParams->isSublist("Contact")){
    Teuchos::ParameterList & contactParams = problemParams->sublist("Contact");
    computeContact = true;
    if(!contactParams.isParameter("Search Radius"))
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Contact parameter \"Search Radius\" not specified.");
    contactSearchRadius = contactParams.get<double>("Search Radius");
    if(!contactParams.isParameter("Search Frequency"))
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Contact parameter \"Search Frequency\" not specified.");
    contactSearchFrequency = contactParams.get<int>("Search Frequency");
  }


}

void Peridigm::Peridigm::initializeMaterials() {

  // Extract problem parameters sublist
  Teuchos::RCP<Teuchos::ParameterList> problemParams = Teuchos::rcp(&(peridigmParams->sublist("Problem")),false);

  // Instantiate material objects
  //! \todo Move creation of material models to material model factory
  TEST_FOR_EXCEPT_MSG(!problemParams->isSublist("Material"), "Material parameters not specified!");
  Teuchos::ParameterList & materialParams = problemParams->sublist("Material");
  Teuchos::ParameterList::ConstIterator it;
  int scalarConstitutiveDataSize = 1; // Epetra barfs if you try to create a vector of zero size
  int vectorConstitutiveDataSize = 1;
  int bondConstitutiveDataSize = 1;
  for(it = materialParams.begin() ; it != materialParams.end() ; it++){
    const string & name = it->first;
    Teuchos::ParameterList & matParams = materialParams.sublist(name);
    Teuchos::RCP<Material> material;
    if(name == "Linear Elastic" || name == "Elastic-Plastic"){
      if(name == "Linear Elastic")
        material = Teuchos::rcp(new LinearElasticIsotropicMaterial(matParams) );
      else if(name == "Elastic-Plastic")
        material = Teuchos::rcp(new IsotropicElasticPlasticMaterial(matParams) );
      materials.push_back( Teuchos::rcp_implicit_cast<Material>(material) );
      // Allocate enough space for the max number of state variables
      if(material->NumScalarConstitutiveVariables() > scalarConstitutiveDataSize)
            scalarConstitutiveDataSize = material->NumScalarConstitutiveVariables();
      if(material->NumVectorConstitutiveVariables() > vectorConstitutiveDataSize)
            vectorConstitutiveDataSize = material->NumVectorConstitutiveVariables();
      if(material->NumBondConstitutiveVariables() > bondConstitutiveDataSize)
            bondConstitutiveDataSize = material->NumBondConstitutiveVariables();
    }
    else {
      string invalidMaterial("Unrecognized material model: ");
      invalidMaterial += name;
      invalidMaterial += ", must be Linear Elastic or Elastic-Plastic";
      TEST_FOR_EXCEPT_MSG(true, invalidMaterial);
    }
  }
  TEST_FOR_EXCEPT_MSG(materials.size() == 0, "No material models created!");

}


