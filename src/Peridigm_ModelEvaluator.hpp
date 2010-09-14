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

#ifndef PERIDIGM_MODELEVALUATOR_HPP
#define PERIDIGM_MODELEVALUATOR_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_LocalMap.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <EpetraExt_ModelEvaluator.h>
#include <Phalanx.hpp>
#include "PHAL_PeridigmTraits.hpp"
#include "Peridigm_AbstractDiscretization.hpp"
#include "materials/Peridigm_Material.hpp"
#include "materials/Peridigm_LinearElasticIsotropicMaterial.hpp"
#include "materials/Peridigm_IsotropicElasticPlasticMaterial.hpp"
#include "contact/Peridigm_ContactModel.hpp"
#include "contact/Peridigm_ShortRangeForceContactModel.hpp"
#include <vector>

namespace PeridigmNS {

  /*! \brief The main ModelEvaluator class; provides the interface 
   * between the driver code and the computational routines.
   */ 
  class ModelEvaluator : public EpetraExt::ModelEvaluator {

  public:

    //! Obsolete constructor
	ModelEvaluator(const Teuchos::RCP<const Epetra_Comm>& comm,
				   const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Constructor
    ModelEvaluator(const Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > materialsList,
                   const Teuchos::RCP<const Epetra_Comm>& comm);

    //! Destructor
	~ModelEvaluator();

    /** \name Overridden from EpetraExt::ModelEvaluator. */
    //@{

    //! Return solution vector map (pure virtual in EpetraExt::ModelEvaluator)
    Teuchos::RCP<const Epetra_Map> get_x_map() const;

    //! Return residual vector map (pure virtual in EpetraExt::ModelEvaluator)
    Teuchos::RCP<const Epetra_Map> get_f_map() const;

    //! Return parameter vector map
    Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;

    //! Return initial x map
    Teuchos::RCP<const Epetra_Vector> get_x_init() const;
    
    //! Return the initial parameter values
    Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;

    //! Return the map for the response function
    Teuchos::RCP<const Epetra_Map> get_g_map(int l) const;

    //! Create InArgs (pure virtual in EpetraExt::ModelEvaluator)
    InArgs createInArgs() const;

    //! Create OutArgs (pure virtual in EpetraExt::ModelEvaluator)
    OutArgs createOutArgs() const;

    //! Evaluate model on InArgs (pure virtual in EpetraExt::ModelEvaluator)
    void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;

    //! Model evaluation that acts directly on the workset
    void evalModel(Teuchos::RCP<PHAL::Workset> workset) const;

    //! Callback function for updating contact.
    virtual void updateContact(Teuchos::RCP<const Epetra_Vector> solverX);

    //! Update the contact neighbor list.  This involves a rebalance.
    virtual void updateContactNeighborList(Teuchos::RCP<const Epetra_Vector> solverX);

    //@}

    //! Update internal history-dependent state information
    void updateState() {};

	//! Return one-dimensional map.
	Teuchos::RCP<const Epetra_Map> getOneDimensionalMap() const;

	//! Return one-dimensional overlap map.
	Teuchos::RCP<const Epetra_Map> getOneDimensionalOverlapMap() const;

    //! Return array of materials
    //! \todo Return reference to array or change to ArrayRCP
    Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > getMaterials() const;

    //! Return array of contact models
    std::vector< Teuchos::RCP<PeridigmNS::ContactModel> > getContactModels() const;

    //! Return neighborlist
    Teuchos::RCP<PeridigmNS::NeighborhoodData> getNeighborhoodData() const;

	//! Return scalar constitutive data
	Teuchos::RCP<const Epetra_MultiVector> getScalarConstitutiveDataOverlap() const;

	//! Return vector constitutive data
	Teuchos::RCP<const Epetra_MultiVector> getVectorConstitutiveDataOverlap() const;

	//! Return bond constitutive data
	Teuchos::RCP<const Epetra_MultiVector> getBondConstitutiveData() const;

  protected:

	void constructEvaluators();

	//! Compute global residual
    void computeGlobalResidual(Teuchos::RCP<const Epetra_Vector>& solverX, 
							   Teuchos::RCP<Epetra_Vector>& solverXDot, 
							   double timeStep) const;

    //! Evaluate response functions
    virtual void 
    evaluateResponses(const Epetra_Vector* x_, 
                      const Epetra_Vector& xdot_,
                      Epetra_Vector& g) const;

	virtual void
	applyBoundaryConditions(const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Epetra map for parameter vector
 	Teuchos::RCP<Epetra_LocalMap> epetraParamMap;

    //! Epetra parameter vector
 	Teuchos::RCP<Epetra_Vector> epetraParamVec;

    //! Supports parameters
    bool supportsP;

    //! Supports response functions
    bool supportsG;

	//! Verbosity flag
	bool verbose;

	// \todo Comments here.
    Teuchos::RCP<const Epetra_Map> oneDimensionalMap;
    Teuchos::RCP<const Epetra_Map> oneDimensionalOverlapMap;
    Teuchos::RCP<const Epetra_Map> threeDimensionalOverlapMap;
    Teuchos::RCP<const Epetra_Map> threeDimensionalTwoEntryMap;
    Teuchos::RCP<const Epetra_Map> secondaryEntryOverlapMap;

	double* bondData;
    Teuchos::RCP<const Epetra_Map> bondMap;

	//! Discretization
	Teuchos::RCP<PeridigmNS::AbstractDiscretization> disc;

    //! Initial positions (vector formatted to interface with the solver)
    Teuchos::RCP<Epetra_Vector> solverInitialX;

    //! Cell volumes
    Teuchos::RCP<Epetra_Vector> cellVolumeOverlap;

    //! Initial positions (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> xOverlap;

    //! Displacement at previously accepted solution (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> uOverlap;

    //! Current displacement (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> yOverlap;

    //! Current velocities (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> vOverlap;

    //! Internal force vector (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> forceOverlap;

    //! Force due to contact (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> contactForceOverlap;

    //! Scalar constitutive data (vector includes ghosted dof)
    Teuchos::RCP<Epetra_MultiVector> scalarConstitutiveDataOverlap;

    //! Bond constitutive data
    Teuchos::RCP<Epetra_MultiVector> bondConstitutiveData;

    //! Vector constitutive data (vector includes ghosted dof)
    Teuchos::RCP<Epetra_MultiVector> vectorConstitutiveDataOverlap;

	/** \brief Importer containing a mapping between the Peridigm overlap
	 *  vectors (xOverlap, vOverlap, etc) and the first entry in a
	 *  vector formatted for the solver (for example, the velocities in 
	 *  the solver's x_dot vector, which contains interleaved velocities
	 *  and accelerations). */
    Teuchos::RCP<Epetra_Import> firstEntryImporter;

	/** \brief Importer containing a mapping between the Peridigm overlap
	 *  vectors (xOverlap, vOverlap, etc) and the second entry in a
	 *  vector formatted for the solver (for example, the accelerations in 
	 *  the solver's x_dot vector, which contains interleaved velocities
	 *  and accelerations). */
	Teuchos::RCP<Epetra_Import> secondEntryImporter;

//     //! Parameter library
//     Teuchos::RCP<ParamLib> paramLib;

	//! Map for combined response functions
    Teuchos::RCP<Epetra_Map> responseMap;

	//! Phalanx field manager
	Teuchos::RCP<PHX::FieldManager<PHAL::PeridigmTraits> > fm;

	//! Output stream
	Teuchos::RCP<Teuchos::FancyOStream> out;

	//! Material models
	//! \todo Use Teuchos::ArrayRCP to store materials?
    Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > materials;

	//! List of neighbors for all locally-owned nodes
	Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData;

	//! Contact models
	Teuchos::RCP<std::vector< Teuchos::RCP<const PeridigmNS::ContactModel> > > contactModels;

	//! List of potential contact neighbors for all locally-owned nodes
	Teuchos::RCP<PeridigmNS::NeighborhoodData> contactNeighborhoodData;

    //! Contact flag
    bool computeContact;

    //! Contact search radius
    double contactSearchRadius;

    //! Contact search frequency
    int contactSearchFrequency;

	//! Number of processors
	int numPID;

	//! Processor ID
	int myPID;

  private:
    
    //! Private to prohibit copying
    ModelEvaluator(const ModelEvaluator&);

    //! Private to prohibit copying
    ModelEvaluator& operator=(const ModelEvaluator&);
  };
}

#endif // PERIDIGM_MODELEVALUATOR_HPP
