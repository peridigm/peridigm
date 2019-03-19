/*! \file Peridigm_BlockBase.hpp */

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

#ifndef PERIDIGM_BLOCKBASE_HPP
#define PERIDIGM_BLOCKBASE_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>

#include <vector>
#include <string>
#include <map>

#include "Peridigm_NeighborhoodData.hpp"
#include "Peridigm_DataManager.hpp"

namespace PeridigmNS {

  class BlockBase {
  public:

    //! Constructor
    BlockBase() : blockName("Undefined"), blockID(-1) {}

    //! Constructor
    BlockBase(std::string blockName_, int blockID_, Teuchos::ParameterList& blockParams_);

    //! Destructor
    ~BlockBase(){}

    //! @name Accessor functions for maps
    //@{
    //! Get the map for scalar data stored at owned points
    Teuchos::RCP<const Epetra_BlockMap> getOwnedScalarPointMap(){ return ownedScalarPointMap; }
    //! Get the map for scalar data stored at owned + overlap points
    Teuchos::RCP<const Epetra_BlockMap> getOverlapScalarPointMap(){ return overlapScalarPointMap; }
    //! Get the map for vector data stored at owned points
    Teuchos::RCP<const Epetra_BlockMap> getOwnedVectorPointMap(){ return ownedVectorPointMap; }
    //! Get the map for vector data stored at owned + overlap points
    Teuchos::RCP<const Epetra_BlockMap> getOverlapVectorPointMap(){ return overlapVectorPointMap; }
    //! Get the map for scalar data stored at owned bonds
    Teuchos::RCP<const Epetra_BlockMap> getOwnedScalarBondMap(){ return ownedScalarBondMap; }
    //@}

    //! Get the neighborhood data
    Teuchos::RCP<PeridigmNS::NeighborhoodData> getNeighborhoodData(){
      return neighborhoodData;
    }

    /*! \brief Initialize the block.
     *
     *  This function will create the block-specific maps, neighborhood list, and DataManager.
     */
    void initialize(Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarPointMap,
                    Teuchos::RCP<const Epetra_BlockMap> globalOverlapScalarPointMap,
                    Teuchos::RCP<const Epetra_BlockMap> globalOwnedVectorPointMap,
                    Teuchos::RCP<const Epetra_BlockMap> globalOverlapVectorPointMap,
                    Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarBondMap,
                    Teuchos::RCP<const Epetra_Vector> globalBlockIds,
                    Teuchos::RCP<const PeridigmNS::NeighborhoodData> globalNeighborhoodData);

    //! Stores a list of field ids that will be added to this block's DataManager.
    void setAuxiliaryFieldIds(std::vector<int> fieldIds){
      auxiliaryFieldIds = fieldIds;
    }

    //! Get the DataManager.
    Teuchos::RCP<PeridigmNS::DataManager> getDataManager(){
      return dataManager;
    }

    //! Get the Block ID
    int getID(){
      return blockID;
    }

    //! Get the block name
    std::string getName(){
      return blockName;
    }

    //! Get the number of points in the block (does not include ghosts)
    int numPoints() {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
        ownedScalarPointMap.is_null(),
        "\n**** Error in BlockBase::numPoints():  Map not set, pointer is null.\n");
      return ownedScalarPointMap->NumMyElements();
    }

    //! Method for accessing data from the DataManager.
    Teuchos::RCP<Epetra_Vector> getData(int fieldId, PeridigmField::Step step){
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
        dataManager.is_null(),
        "\n**** DataManager must be initialized via BlockBase::initializeDataManager() prior to calling BlockBase::getData()\n");
      return dataManager->getData(fieldId, step);
    }

    //! Method for querying the DataManager for the presence of a field spec.
    bool hasData(int fieldId, PeridigmField::Step step) const {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
        dataManager.is_null(),
        "\n**** DataManager must be initialized via BlockBase::initializeDataManager() prior to calling BlockBase::hasData()\n");
      return dataManager->hasData(fieldId, step);
    }

    /*! \brief Import data from the given source vector to the underlying target vector associated with the given field spec.
     *
     *  The intended use case is to import from a non-overlapped vector (i.e., no ghosts) to an overlap vector in
     *  the BlockBase's DataManager.  Note that if the BlockBase does not have space allocated for the given spec and step,
     *  then by design the function is a no-op.
     */
    void importData(Teuchos::RCP<const Epetra_Vector> source, int fieldId, PeridigmField::Step step, Epetra_CombineMode combineMode);

    /*! \brief Export data from the underlying source vector associated with the given field spec to the given target vector.
     *
     *  The intended use case is to export from an overlapped vector (i.e., a vector that includes ghosts) to a non-overlap
     *  (global) vector.  Note that if the BlockBase does not have space allocated for the given spec and step, then by design
     *  the function is a no-op.
     */
    void exportData(Teuchos::RCP<Epetra_Vector> target, int fieldId, PeridigmField::Step step, Epetra_CombineMode combineMode);

    //! Swaps STATE_N and STATE_NP1.
    void updateState(){ dataManager->updateState(); };

    //! Write block data
    void writeBlocktoDisk(std::string blockName, char const * path){ dataManager->writeBlocktoDisk(blockName, path); }

    //! Read block data
    void readBlockfromDisk(std::string blockName, char const * path){ dataManager->readBlockfromDisk(blockName, path); }

  protected:

    /*! \brief Creates the set of block-specific maps.
     *
     *  The block-specific maps are a subset of the global maps.  This function creates the
     *  block-specific maps based on the global maps and a vector containing the block ID
     *  of each element.
     */
    void createMapsFromGlobalMaps(Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarPointMap,
                                  Teuchos::RCP<const Epetra_BlockMap> globalOverlapScalarPointMap,
                                  Teuchos::RCP<const Epetra_BlockMap> globalOwnedVectorPointMap,
                                  Teuchos::RCP<const Epetra_BlockMap> globalOverlapVectorPointMap,
                                  Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarBondMap,
                                  Teuchos::RCP<const Epetra_Vector>   globalBlockIds,
                                  Teuchos::RCP<const PeridigmNS::NeighborhoodData> globalNeighborhoodData);

    //! Create the block-specific neighborhood data.
    Teuchos::RCP<PeridigmNS::NeighborhoodData> createNeighborhoodDataFromGlobalNeighborhoodData(Teuchos::RCP<const Epetra_BlockMap> globalOverlapScalarPointMap,
                                                                                                Teuchos::RCP<const PeridigmNS::NeighborhoodData> globalNeighborhoodData);

    /*! \brief Initialize the data manager.
     *
     *  The DataManager will include all the field specs requested by the material model and
     *  the contact model, as well as those provided by setAuxiliaryFieldIds().
     */
    void initializeDataManager(std::vector<int> fieldIds);

    std::string blockName;
    int blockID;

    //! @name Maps
    //@{
    //! One-dimensional map for owned points.
    Teuchos::RCP<const Epetra_BlockMap> ownedScalarPointMap;
    //! One-dimensional overlap map for owned points and ghosts.
    Teuchos::RCP<const Epetra_BlockMap> overlapScalarPointMap;
    //! Three-dimensional map for owned points.
    Teuchos::RCP<const Epetra_BlockMap> ownedVectorPointMap;
    //! Three-dimensional overlap map for owned points and ghosts.
    Teuchos::RCP<const Epetra_BlockMap> overlapVectorPointMap;
    //! Bond map for owned points; the map contains one element for each owned point, the element length is equal to the number of bonds for that point.
    Teuchos::RCP<const Epetra_BlockMap> ownedScalarBondMap;
    //@}

    //! One-dimensional Importer from global to overlapped vectors
    Teuchos::RCP<const Epetra_Import> oneDimensionalImporter;

    //! One-dimensional Importer from global to overlapped vectors
    Teuchos::RCP<const Epetra_Import> threeDimensionalImporter;

    //! The neighborhood data
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData;

    //! List of auxiliary field specs
    std::vector<int> auxiliaryFieldIds;

    //! The DataManager
    Teuchos::RCP<PeridigmNS::DataManager> dataManager;

    //! The blocks parameterlist sublist
    Teuchos::ParameterList blockParams;
  };

}

#endif // PERIDIGM_BLOCKBASE_HPP
