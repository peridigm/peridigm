/*! \file Peridigm_ModelEvaluator.hpp */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

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
    ModelEvaluator(const Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > materialModels,
                   const Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::ContactModel> > > contactModels,
                   const Teuchos::RCP<const Epetra_Comm>& comm);

    //! Destructor
	virtual ~ModelEvaluator();

    //! Model evaluation that acts directly on the workset
    void evalModel(Teuchos::RCP<PHAL::Workset> workset) const;

    //! Model evaluation that acts directly on the workset
    void evalJacobian(Teuchos::RCP<PHAL::Workset> workset) const;

    //! Update internal history-dependent state information
    void updateState() {};

    //! Return vector of materials
    Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > getMaterialModels() const;

    //! Return vector of contact models
    Teuchos::RCP<std::vector< Teuchos::RCP<const PeridigmNS::ContactModel> > > getContactModels() const;

  protected:

	void constructForceEvaluators();
	void constructJacobianEvaluators();

	//! Phalanx field manager for internal force evaluation
	Teuchos::RCP<PHX::FieldManager<PHAL::PeridigmTraits> > forceFieldManager;

	//! Phalanx field manager for jacobian evaluation
	Teuchos::RCP<PHX::FieldManager<PHAL::PeridigmTraits> > jacobianFieldManager;

	//! Material models
    Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > materialModels;

	//! Contact models
  	Teuchos::RCP<std::vector< Teuchos::RCP<const PeridigmNS::ContactModel> > > contactModels;

    //! Contact flag
    bool analysisHasContact;

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
