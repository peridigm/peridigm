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
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_VerboseObject.hpp>

using namespace std;

//**********************************************************************
template<typename EvalT, typename Traits>
UpdateForceState<EvalT, Traits>::UpdateForceState(const Teuchos::ParameterList& p) :
  parameter_accessible_via_getValue(0.0),
  m_verbose(false),
  m_num_pt(0)
{
  if(p.isParameter("Verbose"))
	 m_verbose = p.get<bool>("Verbose");

  /** \todo Connect parameters to ParamLib if allowing parameters to
   *  be modified by external drivers (e.g. Dakota). */
//   Teuchos::RCP<ParamLib> paramLib = p.get< Teuchos::RCP<ParamLib> >("Parameter Library", Teuchos::null);
//   new PHAL::ParameterEntry<EvalT> ("Coefficient", this, paramLib);

  Teuchos::RCP<PHX::FieldTag> gather_neighborhood_data_field_tag = 
    Teuchos::rcp(new PHX::Tag<ScalarT>("GatherNeighborhoodData", p.get< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout")));

  update_force_state_field_tag = 
    Teuchos::rcp(new PHX::Tag<ScalarT>("UpdateForceState", p.get< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout")));

  this->addEvaluatedField(*update_force_state_field_tag);

  this->setName("UpdateForceState");
}

template<typename EvalT, typename Traits>
void UpdateForceState<EvalT, Traits>::setup_vectors(const Teuchos::ParameterList& p)
{
  /** \todo Is this where we're supposed to access vectors via
   *  the DataLayout?  Currently we're just passing ref-count
   *  pointers to the data directly to the evaluateFields()
   *  method.  Need to research what DemoApps applications
   *  are doing. */
}


//**********************************************************************
template<typename EvalT, typename Traits>
void UpdateForceState<EvalT, Traits>::postRegistrationSetup(
                      typename Traits::SetupData d,
                      PHX::FieldManager<Traits>& fm)
{
}

//**********************************************************************
template<typename EvalT, typename Traits>
void UpdateForceState<EvalT, Traits>::evaluateFields(typename Traits::EvalData cellData)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  if(m_verbose)
	cout << "CHECK inside UpdateForceState::evaluateFields()" << endl;

  const Epetra_Vector& u = *cellData.uOverlap;
  const Epetra_Vector& v = *cellData.vOverlap;
  const double dt = *cellData.timeStep;
  const int numOwnedPoints = cellData.neighborhoodData->NumOwnedPoints();
  const int* ownedIDs = cellData.neighborhoodData->OwnedIDs();
  const int* neighborhoodList = cellData.neighborhoodData->NeighborhoodList();

  PeridigmNS::DataManager& dataManager = *cellData.dataManager;

  // \todo expand bondData to allow for an arbitrary number of bond data per bond
  double* bondData = cellData.bondData.get();
  
  Epetra_Vector& force = *cellData.forceOverlap;

  // handling of material models needs work!
  Teuchos::RCP<const PeridigmNS::Material> material = (*cellData.materials)[0];

  material->updateConstitutiveData(u, 
								   v, 
								   dt, 
								   numOwnedPoints,
								   ownedIDs,
								   neighborhoodList,
								   bondData,
                                   dataManager,
								   force);

  // distribute constitutive data across processors here, if required
  //
  // note:  is not required by any current material model, could potentiall be required
  //        in the future.  optional parallel operations like this could be put in their
  //        own evaluators, and the evaluator could be loaded into the field manager
  //        only if needed.
}
