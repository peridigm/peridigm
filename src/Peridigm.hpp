/*! \file Peridigm.hpp */

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

#ifndef PERIDIGM_HPP
#define PERIDIGM_HPP

#include <vector>

#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include "contact/Peridigm_ContactModel.hpp"
#include "materials/Peridigm_Material.hpp"
#include "Peridigm_AbstractDiscretization.hpp"
#include "Peridigm_ModelEvaluator.hpp"
#include "Peridigm_DataManager.hpp"
#include "Peridigm_OutputManager.hpp"

namespace PeridigmNS {

  class Peridigm {

  public:

    //! Constructor
    Peridigm(const Teuchos::RCP<const Epetra_Comm>& comm,
             const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Instantiate material objects
    void instantiateMaterials();

    //! Initialize material objects
    void initializeMaterials();

    //! Initialize discretization and maps
    void initializeDiscretization();

    //! Apply boundary conditions
    void applyInitialVelocities();

    //! Initialize contact
    void initializeContact();

    //! Initialize the workset
    void initializeWorkset();

    //! Initialize the output manager
    void initializeOutputManager();

    //! Main routine to drive problem solution
    void execute();

    //! Rebalance the mesh
    void rebalance();

    //! Update contact neighborlist; do load rebalance
    void updateContactNeighborList();

    //! Accessor for three-dimensional map
    Teuchos::RCP<const Epetra_BlockMap> getThreeDimensionalMap() { return threeDimensionalMap; }

    //! Destructor
    ~Peridigm(){};

  private:

    //! Parameterlist of entire input deck
    Teuchos::RCP<Teuchos::ParameterList> peridigmParams;

    //! Epetra communicator established by Peridigm_Factory
    Teuchos::RCP<const Epetra_Comm> peridigmComm;

    //! Output stream
    Teuchos::RCP<Teuchos::FancyOStream> out;

    //! Maps for scalar and R3 vector data
    Teuchos::RCP<const Epetra_BlockMap> oneDimensionalMap;
    Teuchos::RCP<const Epetra_BlockMap> oneDimensionalOverlapMap;
    Teuchos::RCP<const Epetra_BlockMap> threeDimensionalMap;
    Teuchos::RCP<const Epetra_BlockMap> threeDimensionalOverlapMap;

    //! Importers and exporters from global to overlapped vectors
    Teuchos::RCP<const Epetra_Import> threeDimensionalMapToThreeDimensionalOverlapMapImporter;
    Teuchos::RCP<const Epetra_Import> oneDimensionalMapToOneDimensionalOverlapMapImporter;

    //! Data structures for bonds
    double* bondData;
    Teuchos::RCP<const Epetra_BlockMap> bondMap;

    //! Contact flag
    bool computeContact;

    //! Contact search radius
    double contactSearchRadius;

    //! Contact search frequency
    int contactSearchFrequency;

    //! Material models
    //! \todo Use Teuchos::ArrayRCP to store materials?
    Teuchos::RCP< std::vector< Teuchos::RCP<const PeridigmNS::Material> > > materials;

    //! Contact models
    Teuchos::RCP< std::vector<Teuchos::RCP<const PeridigmNS::ContactModel> > > contactModels;

    //! Data manager
    Teuchos::RCP<PeridigmNS::DataManager> dataManager;

    //! Global vector for initial positions
    Teuchos::RCP<Epetra_Vector> x;

    //! Global vector for displacement
    Teuchos::RCP<Epetra_Vector> u;

    //! Global vector for velocity
    Teuchos::RCP<Epetra_Vector> v;

    //! Global vector for acceleration
    Teuchos::RCP<Epetra_Vector> a;

    //! Global vector for force
    Teuchos::RCP<Epetra_Vector> force;

    //! Initial positions (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> xOverlap;

    //! Displacement at previously accepted solution (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> uOverlap;

    //! Current velocities (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> vOverlap;

    //! Internal force vector (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> forceOverlap;

    //! Force due to contact (vector includes ghosted dof)
    Teuchos::RCP<Epetra_Vector> contactForceOverlap;

    //! Cell volumes
    Teuchos::RCP<Epetra_Vector> cellVolumeOverlap;

    //! Scalar constitutive data (vector includes ghosted dof)
    int scalarConstitutiveDataSize;
    Teuchos::RCP<Epetra_MultiVector> scalarConstitutiveDataOverlap;

    //! Vector constitutive data (vector includes ghosted dof)
    int vectorConstitutiveDataSize;
    Teuchos::RCP<Epetra_MultiVector> vectorConstitutiveDataOverlap;

    //! Bond constitutive data
    int bondConstitutiveDataSize;
    Teuchos::RCP<Epetra_MultiVector> bondConstitutiveData;

    //! List of neighbors for all locally-owned nodes
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData;

    //! List of potential contact neighbors for all locally-owned nodes
    Teuchos::RCP<PeridigmNS::NeighborhoodData> contactNeighborhoodData;

    //! Workset that is passed to the modelEvaluator
    Teuchos::RCP<PHAL::Workset> workset;

    //! The peridigm model evaluator
    Teuchos::RCP<PeridigmNS::ModelEvaluator> modelEvaluator;

    //! The peridigm output manager
    Teuchos::RCP<PeridigmNS::OutputManager> outputManager;
    //! Description of force state data used by output manager
    Teuchos::RCP<Teuchos::ParameterList> forceStateDesc;

  };
}

#endif // PERIDIGM_HPP
