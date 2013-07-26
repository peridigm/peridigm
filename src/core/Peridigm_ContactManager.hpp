/*! \file Peridigm_ContactManager.hpp */

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

#ifndef PERIDIGM_CONTACTMANAGER_HPP
#define PERIDIGM_CONTACTMANAGER_HPP

#include <vector>
#include <map>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include "Peridigm_ContactBlock.hpp"
#include "contact/Peridigm_ContactModel.hpp"
#include "mesh_input/quick_grid/QuickGridData.h"

// \todo These includes are temporary, remove them.
#include "Peridigm_Discretization.hpp"

namespace PeridigmNS {

/*! \brief Handles contact.
 */
  class ContactManager {
  public:

    // \todo Removed last to parameters to ContactManager constructor.

    //! Constructor.
    ContactManager(const Teuchos::ParameterList& contactParams,
                   Teuchos::RCP<PeridigmNS::Discretization> disc,
                   Teuchos::RCP<Teuchos::ParameterList> peridigmParams);

    //! Initialization routine to allocate and initialize data structures.
    void initialize(Teuchos::RCP<const Epetra_BlockMap> oneDimensionalMap_,
                    Teuchos::RCP<const Epetra_BlockMap> threeDimensionalMap_,
                    Teuchos::RCP<const Epetra_BlockMap> oneDimensionalOverlapMap_,
                    Teuchos::RCP<const Epetra_BlockMap> bondMap_,
                    std::map<std::string, double> blockHorizonValues);

    //! \todo Create map that tracks names of fields in mothership data, use it to simplify data transfers, elliminate the superfluous vector pointers
    
    //! Initial copy of all Peridigm mothership data into contact mothership data.
    void loadAllMothershipData(Teuchos::RCP<Epetra_Vector> blockIds,
                               Teuchos::RCP<Epetra_Vector> volume,
                               Teuchos::RCP<Epetra_Vector> y,
                               Teuchos::RCP<Epetra_Vector> v);

    void loadNeighborhoodData(Teuchos::RCP<PeridigmNS::NeighborhoodData> globalNeighborhoodData);

    void initializeContactBlocks();

    void importData(Teuchos::RCP<Epetra_Vector> volume,
                    Teuchos::RCP<Epetra_Vector> coordinates,
                    Teuchos::RCP<Epetra_Vector> velocity);

    void exportData(Teuchos::RCP<Epetra_Vector> contactForce);

    Teuchos::RCP< std::vector<PeridigmNS::ContactBlock> > getContactBlocks(){ return contactBlocks; };

    void rebalance(int step);

    //! Destructor.
    ~ContactManager(){}

  protected:

    //! Boundary and initial condition parameters
    Teuchos::ParameterList params;

  private:

    //! Compute a parallel decomposion based on the current configuration
    QUICKGRID::Data currentConfigurationDecomp();

    //! Create a rebalanced bond map
    Teuchos::RCP<Epetra_BlockMap> createRebalancedBondMap(Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                          Teuchos::RCP<const Epetra_Import> oneDimensionalMapToRebalancedOneDimensionalMapImporter);

    //! Create a global ID neighbor list in a rebalanced partitioning
    Teuchos::RCP<Epetra_Vector> createRebalancedNeighborGlobalIDList(Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap,
                                                                     Teuchos::RCP<const Epetra_Import> bondMapToRebalancedBondMapImporter);

    //! Create a rebalanced NeighborhoodData object
    Teuchos::RCP<PeridigmNS::NeighborhoodData> createRebalancedNeighborhoodData(Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                                                Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalOverlapMap,
                                                                                Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap,
                                                                                Teuchos::RCP<Epetra_Vector> rebalancedNeighborGlobalIDs);

    //! Fill the contact neighbor information in rebalancedDecomp and populate contactNeighborsGlobalIDs and offProcesorContactIDs
    void contactSearch(Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalMap,
                       Teuchos::RCP<const Epetra_BlockMap> rebalancedBondMap,
                       Teuchos::RCP<const Epetra_Vector> rebalancedNeighborGlobalIDs,
                       QUICKGRID::Data& rebalancedDecomp,
                       Teuchos::RCP< std::map<int, std::vector<int> > > contactNeighborGlobalIDs,
                       Teuchos::RCP< std::set<int> > offProcessorContactIDs);

    //! Create a rebalanced NeighborhoodData object for contact
    Teuchos::RCP<PeridigmNS::NeighborhoodData> createRebalancedContactNeighborhoodData(Teuchos::RCP<std::map<int, std::vector<int> > > contactNeighborGlobalIDs,
                                                                                       Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                                                       Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalOverlapMap);

    //! Global (non-rebalanced) maps used by the Peridigm time integrator
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> bondMap;

    //! Internal maps for contact calculations.
    Teuchos::RCP<const Epetra_BlockMap> oneDimensionalContactMap;
    Teuchos::RCP<const Epetra_BlockMap> threeDimensionalContactMap;
    Teuchos::RCP<const Epetra_BlockMap> oneDimensionalOverlapContactMap;
    Teuchos::RCP<const Epetra_BlockMap> bondContactMap;

    //! Contact mothership multivector for three-dimensional data.
    Teuchos::RCP<Epetra_MultiVector> threeDimensionalContactMothership;

    //! Contact mothership multivector for one-dimensional data.
    Teuchos::RCP<Epetra_MultiVector> oneDimensionalContactMothership;

    //! List of neighbors for all locally-owned nodes, stored in current configuration
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData;

    //! List of neighbors for all locally-owned nodes, stored in current configuration
    Teuchos::RCP<PeridigmNS::NeighborhoodData> contactNeighborhoodData;

    //! Contact search frequency
    int contactRebalanceFrequency;

    //! Contact search radius
    double contactSearchRadius;

    //! Contact models
    std::map< std::string, Teuchos::RCP<const PeridigmNS::ContactModel> > contactModels;

    //! Contact blocks
    Teuchos::RCP< std::vector<PeridigmNS::ContactBlock> > contactBlocks;

    //! Contact block iterator, for convenience
    std::vector<PeridigmNS::ContactBlock>::iterator contactBlockIt;

   //! Global contact vector for block id
    Teuchos::RCP<Epetra_Vector> contactBlockIDs;

    //! Global contact vector for volume
    Teuchos::RCP<Epetra_Vector> contactVolume;

    //! Global contact vector for current position
    Teuchos::RCP<Epetra_Vector> contactY;

    //! Global contact vector for velocity
    Teuchos::RCP<Epetra_Vector> contactV;

    //! Global contact vector for contact force
    Teuchos::RCP<Epetra_Vector> contactContactForce;

    //! Global contact vector for scratch data
    Teuchos::RCP<Epetra_Vector> contactScratch;

    //! Importer for passing one-dmensional data between the mothership vectors and the contact mothership vectors
    Teuchos::RCP<const Epetra_Import> oneDimensionalMothershipToContactMothershipImporter;

    //! Importer for passing three-dmensional data between the mothership vectors and the contact mothership vectors
    Teuchos::RCP<const Epetra_Import> threeDimensionalMothershipToContactMothershipImporter;

    // field ids for all relevant data
    int blockIdFieldId;
    int volumeFieldId;
    int coordinatesFieldId;
    int velocityFieldId;
    int contactForceDensityFieldId;

    // Private to prohibit use.
    ContactManager();

    // Private to prohibit use.
    ContactManager(const ContactManager&);

    // Private to prohibit use.
    ContactManager& operator=(const ContactManager&);
  };
}

#endif // PERIDIGM_CONTACTMANAGER_HPP
