// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#include <Teuchos_TestForException.hpp>
#include <Phalanx_DataLayout.hpp>

//**********************************************************************
template<typename EvalT, typename Traits>
EvaluateForce<EvalT, Traits>::EvaluateForce(Teuchos::ParameterList& p) :
  m_verbose(false),
  m_num_pt(0)
{
  if(p.isParameter("Verbose"))
	 m_verbose = p.get<bool>("Verbose");

  Teuchos::RCP<PHX::FieldTag> update_force_state_data_field_tag = 
	Teuchos::rcp(new PHX::Tag<ScalarT>("UpdateForceState", p.get< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout")));

  this->addDependentField(*update_force_state_data_field_tag);

  evaluate_force_field_tag = 
    Teuchos::rcp(new PHX::Tag<ScalarT>("EvaluateForce",p.get< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout")));

  this->addEvaluatedField(*evaluate_force_field_tag);

  this->setName("EvaluateForce");
}

template<typename EvalT, typename Traits>
void EvaluateForce<EvalT, Traits>::setup_vectors(const Teuchos::ParameterList& p)
{
  // See todo comment in analogous function in PHAL_UpdateForceState_Def.hpp.
}


//**********************************************************************
template<typename EvalT, typename Traits>
void EvaluateForce<EvalT, Traits>::postRegistrationSetup(
                      typename Traits::SetupData d,
                      PHX::FieldManager<Traits>& fm)
{
}

//**********************************************************************
template<typename EvalT, typename Traits>
void EvaluateForce<EvalT, Traits>::evaluateFields(typename Traits::EvalData cellData)
{
  if(m_verbose)
	cout << "CHECK inside EvaluateForce::evaluateFields()\n" << endl;

  const Epetra_Vector& x = *cellData.xOverlap;
  const Epetra_Vector& u = *cellData.uOverlap;
  const Epetra_Vector& v = *cellData.vOverlap;
  const double dt = *cellData.timeStep;
  const Epetra_Vector& cellVolume = *cellData.cellVolumeOverlap;
  const int numOwnedPoints = cellData.neighborhoodData->NumOwnedPoints();
  const int* ownedIDs = cellData.neighborhoodData->OwnedIDs();
  const int* neighborhoodList = cellData.neighborhoodData->NeighborhoodList();

  PeridigmNS::DataManager& dataManager = *cellData.dataManager;

  // \todo expand bondData to allow for an arbitrary number of bond data per bond
  double* bondData = cellData.bondData.get();

  Epetra_MultiVector& vectorConstitutiveData = *cellData.vectorConstitutiveDataOverlap;
  Epetra_Vector& force = *cellData.forceOverlap;
  force.PutScalar(0.0);

  // handling of material models needs work!
  Teuchos::RCP<const PeridigmNS::Material> material = (*cellData.materials)[0];

  material->computeForce(x, 
						 u, 
						 v, 
						 dt, 
						 cellVolume,
						 numOwnedPoints,
						 ownedIDs,
						 neighborhoodList,
						 bondData,
                         dataManager,
						 vectorConstitutiveData,
						 force);
}

