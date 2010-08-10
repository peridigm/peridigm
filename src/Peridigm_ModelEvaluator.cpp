/*! \file Peridigm_ModelEvaluator.hpp */

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

#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_VerboseObject.hpp>
#include "PHAL_Dimension.hpp"
#include "PHAL_FactoryTraits.hpp"
#include "Peridigm_ModelEvaluator.hpp"
#include "Peridigm_DiscretizationFactory.hpp"
#include "PdGridData.h"
#include "PdQuickGrid.h"
#include "PdZoltan.h"
#include <sstream>

Peridigm::ModelEvaluator::ModelEvaluator(const Teuchos::RCP<const Epetra_Comm>& comm,
										 const Teuchos::RCP<Teuchos::ParameterList>& params)
  : supportsP(false),
    supportsG(false),
	verbose(false),
	bondData(0),
    computeContact(false),
    contactSearchRadius(0.0),
    contactSearchFrequency(1),
	numPID(comm->NumProc()),
	myPID(comm->MyPID())
{
  Teuchos::RCP<Teuchos::ParameterList> problemParams = 
	Teuchos::rcp(&(params->sublist("Problem")),false);
  out = Teuchos::VerboseObjectBase::getDefaultOStream();
  verbose = problemParams->get("Verbose", false);

  // Create discretization object
  // The discretization object creates:
  // 1) The epetra maps
  // 2) The vector of initial positions
  Teuchos::RCP<Teuchos::ParameterList> discParams = 
	Teuchos::rcp(&(problemParams->sublist("Discretization")), false);
  Peridigm::DiscretizationFactory discFactory(discParams);
  disc = discFactory.create(comm);

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

  // oneDimensionalMap
  // used for cell volumes and scalar constitutive data
  oneDimensionalMap = disc->getOneDimensionalMap(); //! \todo Is this used?

  // oneDimensionalOverlapMap
  // used for cell volumes and scalar constitutive data
  // includes ghosts
  oneDimensionalOverlapMap = disc->getOneDimensionalOverlapMap();

  // threeDimensionalOverlapMap
  // used for positions, displacements, velocities and vector constitutive data
  // includes ghosts
  threeDimensionalOverlapMap = disc->getThreeDimensionalOverlapMap();
  
  // map for solver unknowns
  // the solver x vector contains current positions and velocities
  // the solver x_dot vector contains velocities and accelerations
  // the order of myGlobalElements is important here, need to allow
  // for import/export with threeDimensionalOverlapMap and 
  // accelerationExportOverlapMap
  // this map must be of type Epetra_Map (not Epetra_BlockMap) so it can be returned by get_x_map()
  threeDimensionalTwoEntryMap = disc->getThreeDimensionalTwoEntryMap();

  // secondaryEntryOverlapMap
  // used for force/acceleration vector and velocity vector
  // allows for import/export into solver's x and x_dot vectors
  secondaryEntryOverlapMap = disc->getSecondaryEntryOverlapMap();

  // bondConstitutiveDataMap
  // a non-overlapping map used for storing constitutive data on bonds
  bondMap = disc->getBondMap();

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
    Teuchos::RCP<Peridigm::Material> material;
	if(name == "Linear Elastic" || name == "Elastic-Plastic"){
      if(name == "Linear Elastic"){
        material = Teuchos::rcp(new LinearElasticIsotropicMaterial(matParams) );
      }
      else if(name == "Elastic-Plastic"){
        material = Teuchos::rcp(new IsotropicElasticPlasticMaterial(matParams) );
      }
	  materials.push_back( Teuchos::rcp_implicit_cast<Peridigm::Material>(material) );
	  // Allocate enough space for the max number of state variables
	  if(material->NumScalarConstitutiveVariables() > scalarConstitutiveDataSize)
		scalarConstitutiveDataSize = material->NumScalarConstitutiveVariables();
	  if(material->NumVectorConstitutiveVariables() > vectorConstitutiveDataSize)
		vectorConstitutiveDataSize = material->NumVectorConstitutiveVariables();
	  if(material->NumBondConstitutiveVariables() > bondConstitutiveDataSize)
		bondConstitutiveDataSize = material->NumBondConstitutiveVariables();
	}
	else{
	  string invalidMaterial("Unrecognized material model: ");
	  invalidMaterial += name;
      invalidMaterial += ", must be Linear Elastic or Elastic-Plastic";
	  TEST_FOR_EXCEPT_MSG(true, invalidMaterial);
	}
  }
  TEST_FOR_EXCEPT_MSG(materials.size() == 0, "No material models created!");

  // Instantiate contact models
  //! \todo Move creation of contact models to contact model factory
  if(computeContact){
    Teuchos::ParameterList & contactParams = problemParams->sublist("Contact");
    Teuchos::ParameterList::ConstIterator it;
    for(it = contactParams.begin() ; it != contactParams.end() ; it++){
      const string & name = it->first;
      if(contactParams.isSublist(name)){
        Teuchos::ParameterList & contactModelParams = contactParams.sublist(name);
        // Add the horizon to the contact model parameters, if needed
        if(!contactModelParams.isParameter("Horizon"))
          contactModelParams.set("Horizon", disc->getHorizon());
        Teuchos::RCP<Peridigm::ContactModel> contactModel;
        if(name == "Short Range Force"){
          contactModel = Teuchos::rcp(new ShortRangeForceContactModel(contactModelParams) );
          contactModels.push_back( Teuchos::rcp_implicit_cast<Peridigm::ContactModel>(contactModel) );
        }
        else{
          string invalidContactModel("Unrecognized contact model: ");
          invalidContactModel += name;
          invalidContactModel += ", must be Short Range Force";
          TEST_FOR_EXCEPT_MSG(true, invalidContactModel);
        }
      }
    }
  }

  // create the importer and exporter objects
  firstEntryImporter = Teuchos::rcp(new Epetra_Import(*threeDimensionalOverlapMap, *threeDimensionalTwoEntryMap));
  secondEntryImporter = Teuchos::rcp(new Epetra_Import(*secondaryEntryOverlapMap, *threeDimensionalTwoEntryMap));

  // Allocate epetra vectors
  solverInitialX = disc->getSolverInitialX();

  // Apply initial velocity boundary conditions
  if(problemParams->isSublist("Boundary Conditions")){
	Teuchos::RCP<Teuchos::ParameterList> boundaryConditionParams = 
	  Teuchos::rcp(&(problemParams->sublist("Boundary Conditions")), false);
	applyBoundaryConditions(boundaryConditionParams);
  }

  // Create x, u, and v vectors
  xOverlap = Teuchos::rcp(new Epetra_Vector(*threeDimensionalOverlapMap));
  xOverlap->Import(*solverInitialX, *firstEntryImporter, Insert);
  uOverlap = Teuchos::rcp(new Epetra_Vector(*threeDimensionalOverlapMap));
  vOverlap = Teuchos::rcp(new Epetra_Vector(*secondaryEntryOverlapMap));
  yOverlap = Teuchos::rcp(new Epetra_Vector(*threeDimensionalOverlapMap));

  // Get the cell volumes and put them in the cellVolumeOverlap vector
  cellVolumeOverlap = Teuchos::rcp(new Epetra_Vector(*oneDimensionalOverlapMap));
  Epetra_Import oneDimensionalMapToOneDimensionalOverlapMapImporter(*oneDimensionalOverlapMap, *oneDimensionalMap);
  cellVolumeOverlap->Import(*(disc->getCellVolume()), oneDimensionalMapToOneDimensionalOverlapMapImporter, Insert);

  // containers for constitutive data
  scalarConstitutiveDataOverlap = Teuchos::rcp(new Epetra_MultiVector(*oneDimensionalOverlapMap, scalarConstitutiveDataSize));
  scalarConstitutiveDataOverlap->PutScalar(0.0);
  vectorConstitutiveDataOverlap = Teuchos::rcp(new Epetra_MultiVector(*threeDimensionalOverlapMap, vectorConstitutiveDataSize));
  vectorConstitutiveDataOverlap->PutScalar(0.0);
  bondConstitutiveData = Teuchos::rcp(new Epetra_MultiVector(*bondMap, bondConstitutiveDataSize));
  bondConstitutiveData->PutScalar(0.0);

  // container for accelerations
  forceOverlap = Teuchos::rcp(new Epetra_Vector(*secondaryEntryOverlapMap));  

  // get the neighborlist from the discretization
  neighborhoodData = disc->getNeighborhoodData();

  // container for accelerations due to contact
  if(computeContact){
    contactForceOverlap = Teuchos::rcp(new Epetra_Vector(*secondaryEntryOverlapMap));  
    contactNeighborhoodData = Teuchos::rcp(new Peridigm::NeighborhoodData);
    updateContactNeighborList(solverInitialX);
  }

  // container for bond damage
  bondData = new double[disc->getNumBonds()];
  for(unsigned int i=0; i<disc->getNumBonds(); i++)
	bondData[i] = 0.0;

  // Initialize material models
  std::vector< Teuchos::RCP<Peridigm::Material> >::const_iterator matIt;
  for(matIt = materials.begin() ; matIt != materials.end() ; matIt++){
    double dt = 0.0;
    (*matIt)->initialize(*xOverlap,
                         *uOverlap,
                         *vOverlap,
                         dt,
                         *cellVolumeOverlap,
                         neighborhoodData->NumOwnedPoints(),
                         neighborhoodData->OwnedIDs(),
                         neighborhoodData->NeighborhoodList(),
                         bondData,
                         *scalarConstitutiveDataOverlap,
                         *vectorConstitutiveDataOverlap,
                         *bondConstitutiveData,
                         *forceOverlap);
  }

  // construct evaluators
  constructEvaluators(problemParams);
  fm->postRegistrationSetup(NULL);

  /** \todo If allowing access to member data (e.g. allowing
   *       Dakota to solve optimization problems), then 
   *       create hooks to a ParamLib here.
   */
}

Peridigm::ModelEvaluator::~ModelEvaluator(){
  if(bondData != 0)
	delete[] bondData;
}

void 
Peridigm::ModelEvaluator::constructEvaluators(const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using std::vector;
  using std::map;
  using PHX::DataLayout;
  using PHX::MDALayout;

  // Create a list of evaluators to build

  // Each field evaluator will be associated with one or more data
  // types and data layouts.

  map<string, RCP<ParameterList> > evaluatorsToBuild;

  // Create a dummy data layout that can be associated with an evaluator if needed
  //! \todo Determine how the field manager's data layout scheme works and set it up properly.
  RCP<DataLayout> dummy = rcp(new MDALayout<Dummy>(0));

  { // Update force state
  RCP<ParameterList> p = rcp(new ParameterList);
  int type = FactoryTraits<PHAL::PeridigmTraits>::id_update_force_state;
  p->set<int>("Type", type); 
  p->set<bool>("Verbose", verbose);

  // Associate the ParameterLibrary with the evaluator
  // this grants outside access (e.g. Dakota) to the parameters in the evaluator
//   p->set<RCP<ParamLib> >("Parameter Library", paramLib);

  // set the coefficient to the value in the parameter list
  // the parameter list could have been read from file or set to a default value
//   p->set<double>("Coefficient", params->get("Coefficient",1.0)); 
//   p->set<Teuchos::ParameterList*>("Parameter List", &paramList);

  p->set< RCP<DataLayout> >("Dummy Data Layout", dummy);

  evaluatorsToBuild["UpdateForceState"] = p;
  }

  { // Evaluate force
  RCP<ParameterList> p = rcp(new ParameterList);
  int type = FactoryTraits<PHAL::PeridigmTraits>::id_evaluate_force;
  p->set<int>("Type", type); 
  p->set<bool>("Verbose", verbose);

  p->set< RCP<DataLayout> >("Dummy Data Layout", dummy);

  evaluatorsToBuild["EvaluateForce"] = p;
  }

  // Contact
  if(computeContact){
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = FactoryTraits<PHAL::PeridigmTraits>::id_contact;
    p->set<int>("Type", type); 
    p->set<bool>("Verbose", verbose);

    p->set< RCP<DataLayout> >("Dummy Data Layout", dummy);

    evaluatorsToBuild["Contact"] = p;
  }

  // Build Field Evaluators for each evaluation type
  PHX::EvaluatorFactory<PHAL::PeridigmTraits, FactoryTraits<PHAL::PeridigmTraits> > factory;
  RCP< vector< RCP<PHX::Evaluator_TemplateManager<PHAL::PeridigmTraits> > > > evaluators;
  evaluators = factory.buildEvaluators(evaluatorsToBuild);

  // Create a FieldManager
  fm = Teuchos::rcp(new PHX::FieldManager<PHAL::PeridigmTraits>);

  // Register all Evaluators with the field manager
  PHX::registerEvaluators(evaluators, *fm);

  // List EvaluateForce as a required field for the field manager (other evaluators will be
  // called as perscribed by dependencies).
  PHX::Tag<PHAL::PeridigmTraits::Residual::ScalarT> evaluate_force_tag("EvaluateForce", dummy);
  fm->requireField<PHAL::PeridigmTraits::Residual>(evaluate_force_tag);

  // Require the contact force evaluation
  if(computeContact){
    PHX::Tag<PHAL::PeridigmTraits::Residual::ScalarT> contact_tag("Contact", dummy);
    fm->requireField<PHAL::PeridigmTraits::Residual>(contact_tag);
  }
}

Teuchos::RCP<const Epetra_Map>
Peridigm::ModelEvaluator::get_x_map() const
{
  return threeDimensionalTwoEntryMap;
}

Teuchos::RCP<const Epetra_Map>
Peridigm::ModelEvaluator::get_f_map() const
{
  return threeDimensionalTwoEntryMap;
}

Teuchos::RCP<const Epetra_Map>
Peridigm::ModelEvaluator::get_p_map(int l) const
{
//  TEST_FOR_EXCEPTION(supportsP == false, 
//                     Teuchos::Exceptions::InvalidParameter,
//                     std::endl << 
//                     "Error!  Peridigm::ModelEvaluator::get_p_map():  " <<
//                     "No parameters have been supplied.  " <<
//                     "Supplied index l = " << l << std::endl);
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  Peridigm::ModelEvaluator::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " << 
                     l << std::endl);

  return epetraParamMap;
}

Teuchos::RCP<const Epetra_Map>
Peridigm::ModelEvaluator::get_g_map(int l) const
{
  TEST_FOR_EXCEPTION(supportsG == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  Peridigm::ModelEvaluator::get_g_map():  " <<
                     "No response functions have been supplied.  " <<
                     "Supplied index l = " << l << std::endl);
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  Peridigm::ModelEvaluator::get_g_map() only " <<
                     " supports 1 response vector.  Supplied index l = " << 
                     l << std::endl);

  return responseMap;
}

Teuchos::RCP<const Epetra_Vector>
Peridigm::ModelEvaluator::get_x_init() const
{
  return solverInitialX;
}

Teuchos::RCP<const Epetra_Vector>
Peridigm::ModelEvaluator::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  Peridigm::ModelEvaluator::get_p_init() only " <<
                     " supports 1 parameter vector.  Supplied index l = " << 
                     l << std::endl);
  
  return epetraParamVec;
}

Teuchos::RCP<Peridigm::NeighborhoodData>
Peridigm::ModelEvaluator::getNeighborhoodData() const
{
  return neighborhoodData;
}

Teuchos::RCP<const Epetra_MultiVector>
Peridigm::ModelEvaluator::getScalarConstitutiveDataOverlap() const
{
  return scalarConstitutiveDataOverlap;
}

Teuchos::RCP<const Epetra_MultiVector>
Peridigm::ModelEvaluator::getVectorConstitutiveDataOverlap() const
{
  return vectorConstitutiveDataOverlap;
}

Teuchos::RCP<const Epetra_MultiVector>
Peridigm::ModelEvaluator::getBondConstitutiveData() const
{
  return bondConstitutiveData;
}

Teuchos::RCP<const Epetra_Map>
Peridigm::ModelEvaluator::getOneDimensionalMap() const
{
  return oneDimensionalMap;
}

Teuchos::RCP<const Epetra_Map>
Peridigm::ModelEvaluator::getOneDimensionalOverlapMap() const
{
  return oneDimensionalOverlapMap;
}

std::vector< Teuchos::RCP<Peridigm::Material> >
Peridigm::ModelEvaluator::getMaterials() const
{
  return materials;
}

EpetraExt::ModelEvaluator::InArgs
Peridigm::ModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x, true);
  inArgs.setSupports(IN_ARG_t, true);

  if (supportsP)
    inArgs.set_Np(1); // 1 vector of parameters of length of epetraParamVec
  else
    inArgs.set_Np(0);

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Peridigm::ModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(OUT_ARG_f, true);
  
  if (supportsP && supportsG) {
    outArgs.set_Np_Ng(1, 1); // 1 parameter vector, 1 response vector
	//! \todo Resolve issue related to Dakota expecting DgDp support if g and p are defined.
	//MLP: DgDp is requested by Dakota. Thyra is (apparently) assuming that if we have g and p
	//     then we must have DgDp, and throws error if we do not support DgDp
    outArgs.setSupports(OUT_ARG_DgDp, 0, 0, 
                        DerivativeSupport(DERIV_MV_BY_COL));

	//! \todo Resolve issue related to Dakota expecting DgDx support if g and x are defined.
	//MLP: Ditto with DgDx
    outArgs.setSupports(OUT_ARG_DgDx, 0, DerivativeSupport(DERIV_MV_BY_COL));
  }
  else if (supportsP)
    outArgs.set_Np_Ng(1, 0); // 1 parameter vector
  else if (supportsG)
    outArgs.set_Np_Ng(0, 1); // 1 response vector
  else
    outArgs.set_Np_Ng(0, 0);

  //! \todo Resolve issue related to setting for DfDp support.
  // THIS SHOULD BE COMMENTED OUT; STATEMENT HERE TO AVOID THYRA BREAKING
  // REPORT BUG TO ROSS 
//   outArgs.setSupports(OUT_ARG_DfDp, 0, DerivativeSupport(DERIV_MV_BY_COL));

  return outArgs;
}

void 
Peridigm::ModelEvaluator::evalModel(const InArgs& inArgs, 
									const OutArgs& outArgs) const
{
  // Note:  This is the call by which Rythmos 
  //        interacts with the ModelEvaluator.
  //        We are only handling one case here,
  //        the case where Rythmos gives the
  //        ModelEvaluator x and asks for x_dot.
  //        In general, there are many other cases,
  //        for example, in a quasistatics
  //        analysis or when performing optimization
  //        or UQ.  See the ModelEvaluators in 
  //        DemoApps for examples of this.

  // In our case:
  //   inArgs.get_x()   returns the unknown vector, 
  //                    which contains both the positions
  //                    x and their velocities x_dot.
  //   outArgs.get_f()  returns the time derivative of
  //                    what we're storing in x.
  // This is confusing because we're solving
  // a second-order PDE as a system of two
  // first-order PDEs.

  // In other cases:
  //   We'd call the appropriate function, e.g.
  //   computeTangent, based on the content of
  //   InArgs and OutArgs.
  //   There's currently no point in checking
  //   the content of InArgs and OutArgs
  //   because we know we're solving only one
  //   type of problem.

  //  y = x + u + v*dt
  //  u = y - x - v*dt

  //! \todo Fix handling of time step, this is awful!
  static double previousTime = 0.0;
  static double currentTime = 0.0;
  previousTime = currentTime;
  currentTime = inArgs.get_t();
  double timeStep = currentTime - previousTime;

  Teuchos::RCP<const Epetra_Vector> solverX = inArgs.get_x();
  Teuchos::RCP<Epetra_Vector> solverXDot = outArgs.get_f();

  yOverlap->Import(*solverX, *firstEntryImporter, Insert);
  vOverlap->Import(*solverX, *secondEntryImporter, Insert);

  //! \todo Use Epetra vector operations here, if possible, or possibly BLAS.
  for(int i=0 ; i<uOverlap->MyLength() ; ++i){
	(*uOverlap)[i] = (*yOverlap)[i] - (*xOverlap)[i] - (*vOverlap)[i]*timeStep;
  }

  computeGlobalResidual(solverX, solverXDot, timeStep);
}

void
Peridigm::ModelEvaluator::computeGlobalResidual(Teuchos::RCP<const Epetra_Vector>& solverX, 
												Teuchos::RCP<Epetra_Vector>& solverXDot, 
												double timeStep) const
{ 
  // set up workset
  PHAL::Workset workset;
  workset.xOverlap = xOverlap;
  workset.uOverlap = uOverlap;
  workset.vOverlap = vOverlap;
  workset.forceOverlap = forceOverlap;
  workset.contactForceOverlap = contactForceOverlap;
  workset.timeStep = Teuchos::RCP<double>(&timeStep, false);
  workset.cellVolumeOverlap = cellVolumeOverlap;
  workset.neighborhoodData = neighborhoodData;
  workset.contactNeighborhoodData = contactNeighborhoodData;
  workset.bondData = Teuchos::RCP<double>(bondData, false);
  workset.scalarConstitutiveDataOverlap = scalarConstitutiveDataOverlap;
  workset.vectorConstitutiveDataOverlap = vectorConstitutiveDataOverlap;
  workset.bondConstitutiveData = bondConstitutiveData;
  workset.materials = materials;
  workset.contactModels = contactModels;
  workset.myPID = myPID;

  // call field manager with workset
  fm->evaluateFields<PHAL::PeridigmTraits::Residual>(workset);

  // add the internal force and the contact force
  if(computeContact){
    forceOverlap->Update(1.0, *contactForceOverlap, 1.0);
  }

  // convert force densities to accelerations
  // \todo Expand to handle multiple materials (would be nice if we could scale the blocks independently)
  double density = materials[0]->Density();
  forceOverlap->Scale(1.0/density);

  // scatter add the forces into the solver's container
  solverXDot->Export(*forceOverlap, *secondEntryImporter, Add);

  // copy the velocities into solverXDot
  // note: don't think we can use scatter here because 
  // underlying maps are inconsistent
  for(int i=0 ; i<oneDimensionalMap->NumMyElements() ; ++i){
	(*solverXDot)[i*6+0] = (*solverX)[i*6+3];
	(*solverXDot)[i*6+1] = (*solverX)[i*6+4];
	(*solverXDot)[i*6+2] = (*solverX)[i*6+5];
  }
}

void
Peridigm::ModelEvaluator::evaluateResponses(const Epetra_Vector* xdot,
											const Epetra_Vector& x,
											Epetra_Vector& g) const
{
  // This function defines the respose, as in the case of
  // an optimization problem where we want to minimize a
  // respose function.  For example, Dakota could interface
  // with this ModelEvaluator and minmize this response
  // function by altering parameters
  TEST_FOR_EXCEPT(false);
}

void
Peridigm::ModelEvaluator::updateContact(Teuchos::RCP<const Epetra_Vector> solverX)
{
  if(!computeContact)
    return;

  static int step = 0;
  step += 1;

  if(step%contactSearchFrequency == 0)
    updateContactNeighborList(solverX);
}

void
Peridigm::ModelEvaluator::updateContactNeighborList(Teuchos::RCP<const Epetra_Vector> solverX)
{
  // initial implementation works in serial only
  TEST_FOR_EXCEPT_MSG(numPID != 1, "Contact is currently not enabled in parallel.\n");

  // Create a decomp object and fill necessary data for rebalance
  int myNumElements = oneDimensionalMap->NumMyElements();
  int dimension = 3;
  PdGridData decomp = PdQuickGrid::allocatePdGridData(myNumElements, dimension);

  // fill myGlobalIDs
  shared_ptr<int> myGlobalIDs(new int[myNumElements], PdQuickGrid::Deleter<int>());
  int* myGlobalIDsPtr = myGlobalIDs.get();
  int* gIDs = oneDimensionalMap->MyGlobalElements();
  for(int i=0 ; i<myNumElements ; ++i){
    myGlobalIDsPtr[i] = gIDs[i];
  }
  decomp.myGlobalIDs = myGlobalIDs;
  
  // fill myX and cellVolume
  shared_ptr<double> myX(new double[myNumElements*dimension], PdQuickGrid::Deleter<double>());
  double* myXPtr = myX.get();
  double* solverXPtr;
  solverX->ExtractView(&solverXPtr);
  shared_ptr<double> cellVolume(new double[myNumElements], PdQuickGrid::Deleter<double>());
  double* cellVolumePtr = cellVolume.get();
  double* cellVolumeOverlapPtr;
  cellVolumeOverlap->ExtractView(&cellVolumeOverlapPtr);
  for(int i=0 ; i<myNumElements ; ++i){
    int oneDimensionalMapGlobalID = myGlobalIDsPtr[i];
    int oneDimensionalOverlapMapLocalID = oneDimensionalOverlapMap->LID(oneDimensionalMapGlobalID);
    int threeDimensionalTwoEntryMapGlobalID = oneDimensionalMapGlobalID*3;
    int threeDimensionalTwoEntryMapLocalID = threeDimensionalTwoEntryMap->LID(threeDimensionalTwoEntryMapGlobalID);
    myXPtr[i*3] = solverXPtr[threeDimensionalTwoEntryMapLocalID];
    myXPtr[i*3+1] = solverXPtr[threeDimensionalTwoEntryMapLocalID+1];
    myXPtr[i*3+2] = solverXPtr[threeDimensionalTwoEntryMapLocalID+2];
    cellVolumePtr[i] = cellVolumeOverlapPtr[oneDimensionalOverlapMapLocalID];
  }  
  decomp.myX = myX;
  decomp.cellVolume = cellVolume;

  // rebalance
  decomp = getLoadBalancedDiscretization(decomp);

  // big todo: shuffle data around based on decomp

  // execute contact search
  decomp = createAndAddNeighborhood(decomp, contactSearchRadius);

  // Copy the data in decomp into the contact neighbor list
  // Do not include points that are bonded
  
  vector<int> contactOwnedIDs;
  vector<int> contactNeighborhoodPtr;
  vector<int> contactNeighborhoodList;

  int searchListIndex = 0;
  int searchNumPoints = decomp.numPoints;
  int* searchNeighborhood = decomp.neighborhood.get();

  for(int iLID=0 ; iLID<searchNumPoints ; ++iLID){

    // find the cells that are bonded to this cell
    // store the corresponding local IDs in bondedNeighbors, which is a stl::list
    int* bondedNeighborhoodList = neighborhoodData->NeighborhoodList();
    int bondedListIndex = neighborhoodData->NeighborhoodPtr()[iLID];
    int numBondedNeighbors = bondedNeighborhoodList[bondedListIndex++];
    list<int> bondedNeighbors; // \todo reserve space here
    for(int i=0 ; i<numBondedNeighbors ; ++i){
      bondedNeighbors.push_back(bondedNeighborhoodList[bondedListIndex++]);
    }

    // loop over the cells found by the contact search
    // retain only those cells that are not bonded
	int searchNumNeighbors = searchNeighborhood[searchListIndex++];

    list<int>::iterator it;
    bool hasContact = false;
    int currentContactNeighborhoodPtr = 0;
	for(int iNeighbor=0 ; iNeighbor<searchNumNeighbors ; ++iNeighbor){
	  int localNeighborID = searchNeighborhood[searchListIndex++];
      it = find(bondedNeighbors.begin(), bondedNeighbors.end(), localNeighborID);
      if(it == bondedNeighbors.end()){
        if(!hasContact){
          hasContact = true;
          contactOwnedIDs.push_back(iLID);
          currentContactNeighborhoodPtr = contactNeighborhoodList.size();
          contactNeighborhoodPtr.push_back(currentContactNeighborhoodPtr);
          contactNeighborhoodList.push_back(1);
        }
        else{
          contactNeighborhoodList[currentContactNeighborhoodPtr] += 1;
        }
        contactNeighborhoodList.push_back(localNeighborID);
      }
	}
  }

  TEST_FOR_EXCEPT_MSG(contactNeighborhoodPtr.size() != contactOwnedIDs.size(),
                      "Error, contactOwnedIDs and contactNeighborhoodPtr are different sizes in ModelEvaluator::updateContactNeighborList().\n");

  // copy the contact neighbor data into contactNeighborData
  contactNeighborhoodData->SetNumOwned(contactOwnedIDs.size());
  memcpy(contactNeighborhoodData->OwnedIDs(), 
		 &contactOwnedIDs[0],
		 contactOwnedIDs.size()*sizeof(int));
  memcpy(contactNeighborhoodData->NeighborhoodPtr(), 
		 &contactNeighborhoodPtr[0],
		 contactOwnedIDs.size()*sizeof(int));
  contactNeighborhoodData->SetNeighborhoodListSize(contactNeighborhoodList.size());
  memcpy(contactNeighborhoodData->NeighborhoodList(),
		 &contactNeighborhoodList[0],
		 contactNeighborhoodList.size()*sizeof(int));
}

void
Peridigm::ModelEvaluator::applyBoundaryConditions(const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  //! \todo Need a better way to parse boundary conditions.
  //! \todo Add assertions to boundary contition parsing

  Teuchos::ParameterList::ConstIterator it;

  // get the node sets
  map< string, vector<int> > nodeSets;
  for(it = params->begin() ; it != params->end() ; it++){
	const string & name = it->first;
	size_t position = name.find("Node Set");
	if(position != string::npos){
	  stringstream ss(Teuchos::getValue<string>(it->second));
	  vector<int> nodeList;
	  int nodeID;
	  while(ss.good()){
		ss >> nodeID;
		nodeList.push_back(nodeID);
	  }
	  nodeSets[name] = nodeList;
	}
  }

  // apply the initial conditions
  for(it = params->begin() ; it != params->end() ; it++){
	const string & name = it->first;
	size_t position = name.find("Initial Velocity");
	if(position != string::npos){
	  Teuchos::ParameterList & bcParams = Teuchos::getValue<Teuchos::ParameterList>(it->second);
	  string nodeSet = bcParams.get<string>("Node Set");
	  string type = bcParams.get<string>("Type");
	  string coordinate = bcParams.get<string>("Coordinate");
	  double value = bcParams.get<double>("Value");

	  int coord = 0;
	  if(coordinate == "y" || coordinate == "Y")
		coord = 1;
	  if(coordinate == "z" || coordinate == "Z")
		coord = 2;

	  // apply initial velocity boundary conditions
	  // to locally-owned nodes
	  vector<int> & nodeList = nodeSets[nodeSet];
	  for(unsigned int i=0 ; i<nodeList.size() ; i++){
		int localNodeID = oneDimensionalMap->LID(nodeList[i]);
		if(localNodeID != -1){
		  // solverInitialX vector is formatted for interface with the solver
		  // the format is [X1, Y1, Z1, Vx1, Vy1, Vz1, X2, Y2, Z2, Vx2, Vy2, Vz2, ...]
		  int index = localNodeID*6 + 3 + coord;
		  (*solverInitialX)[index] = value;
		}
	  }
	}
  }
}
