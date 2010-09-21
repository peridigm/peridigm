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
  class ModelEvaluator {

  public:

    //! Constructor
    ModelEvaluator(const Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > materialsList,
                   const Teuchos::RCP<const Epetra_Comm>& comm);

    //! Destructor
	virtual ~ModelEvaluator();

    //! Model evaluation that acts directly on the workset
    void evalModel(Teuchos::RCP<PHAL::Workset> workset) const;

    //! Callback function for updating contact.
//     virtual void updateContact(Teuchos::RCP<const Epetra_Vector> solverX);

    //! Update the contact neighbor list.  This involves a rebalance.
//     virtual void updateContactNeighborList(Teuchos::RCP<const Epetra_Vector> solverX);

    //@}

    //! Update internal history-dependent state information
    void updateState() {};

    //! Return vector of materials
    Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > getMaterials() const;

    //! Return vector of contact models
//     std::vector< Teuchos::RCP<PeridigmNS::ContactModel> > getContactModels() const;

  protected:

	void constructEvaluators();

	//! Phalanx field manager
	Teuchos::RCP<PHX::FieldManager<PHAL::PeridigmTraits> > fm;

	//! Material models
    Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > materials;

	//! Contact models
//  	Teuchos::RCP<std::vector< Teuchos::RCP<const PeridigmNS::ContactModel> > > contactModels;

	//! List of potential contact neighbors for all locally-owned nodes
// 	Teuchos::RCP<PeridigmNS::NeighborhoodData> contactNeighborhoodData;

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

    //! Verbosity flag
    bool verbose;

  private:
    
    //! Private to prohibit copying
    ModelEvaluator(const ModelEvaluator&);

    //! Private to prohibit copying
    ModelEvaluator& operator=(const ModelEvaluator&);
  };
}

#endif // PERIDIGM_MODELEVALUATOR_HPP
