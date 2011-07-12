/*! \file Peridigm_Block.hpp */

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

#ifndef PERIDIGM_BLOCK_HPP
#define PERIDIGM_BLOCK_HPP

/*

Items that correspond to particular blocks:

MaterialModel
ContactModel
DataManager
Maps:       ownedScalarPointMap, overlapScalarPointMap, ownedVectorPointMap, overlapVectorPointMap, ownedScalarBondMap
Importers:  oneDimensionalMapToOneDimensionalOverlapMapImporter, threeDimensionalMapToThreeDimensionalOverlapMapImporter;

FieldSpecs?

Is there any reason for a global overlap map?  I don't think so.
The global maps are just the globalOneDimensionalMap and the globalThreeDimensionalMap.
I'm not even sure there is a need for a globalOneDimensionalMap.  The mothership
Epetra_MultiVector is built with a single map, the globalOneDimensionalMap.

Need to store volume as mothership vector, because it is needed for rebalance?

Code flow:
1)  Have discretizations supply maps and importers.  The discretizations will also return the global 1D and 3D maps.  Peridigm will never know about the block-specific maps.
2)  Peridigm must create and set the material model and contact model.  Block should have setMaterialModel() and setContactModel() methods.
3)  After the specs are known, Peridigm will create the dataManager (but the block will already have the maps, so this will be pretty clean).  Block should have setSpecs() methods, will already have access to material model and contact model specs, needs only global specs.  Perhaps and addSpec() method and finalizeDataManager() methods.
4)  Block should have spec API similar to current DataManager API so that the output manager can work directly on the blocks.

*/

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>

#include <string>

#include "Peridigm_NeighborhoodData.hpp"
#include "Peridigm_DataManager.hpp"
#include "contact/Peridigm_ContactModel.hpp"
#include "materials/Peridigm_Material.hpp"

namespace PeridigmNS {

  class Block {
  public:

    //! Constructor
    Block() : blockName("Undefined") {}

    //! Constructor
    Block(std::string blockName_) : blockName(blockName_) {}

    //! Destructor
    ~Block(){}

    //! @name Accessor functions for maps
    //@{
    //! Get the map for scalar data stored at owned points
    Teuchos::RCP<const Epetra_BlockMap> getOwnedScalarPointMap(){ return ownedScalarPointMap; }
    //! Set the map for scalar data stored at owned points
    void setOwnedScalarPointMap(Teuchos::RCP<const Epetra_BlockMap> map){ ownedScalarPointMap = map; }
    //! Get the map for scalar data stored at owned + overlap points
    Teuchos::RCP<const Epetra_BlockMap> getOverlapScalarPointMap(){ return overlapScalarPointMap; }
    //! Set the map for scalar data stored at owned + overlap points
    void setOverlapScalarPointMap(Teuchos::RCP<const Epetra_BlockMap> map){ overlapScalarPointMap = map; }
    //! Get the map for vector data stored at owned points
    Teuchos::RCP<const Epetra_BlockMap> getOwnedVectorPointMap(){ return ownedVectorPointMap; }
    //! Set the map for vector data stored at owned points
    void setOwnedVectorPointMap(Teuchos::RCP<const Epetra_BlockMap> map){ ownedVectorPointMap = map; }
    //! Get the map for vector data stored at owned + overlap points
    Teuchos::RCP<const Epetra_BlockMap> getOverlapVectorPointMap(){ return overlapVectorPointMap; }
    //! Set the map for vector data stored at owned + overlap points
    void setOverlapVectorPointMap(Teuchos::RCP<const Epetra_BlockMap> map){ overlapVectorPointMap = map; }
    //! Get the map for scalar data stored at owned bonds
    Teuchos::RCP<const Epetra_BlockMap> getOwnedScalarBondMap(){ return ownedScalarBondMap; }
    //! Set the map for scalar data stored at owned bonds
    void setOwnedScalarBondMap(Teuchos::RCP<const Epetra_BlockMap> map){ ownedScalarBondMap = map; }
    //@}

    //! Set the neighborhood data
    void setNeighborhoodData(Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData_){ neighborhoodData = neighborhoodData_; }

    //! Set the material model
    void setMaterialModel(Teuchos::RCP<const PeridigmNS::Material> materialModel_){ materialModel = materialModel_; }

    //! Set the contact model
    void setContactModel(Teuchos::RCP<const PeridigmNS::ContactModel> contactModel_){ contactModel = contactModel_; }

    /*! \brief Initialize the data manager.
     *
     *  The DataManager will include all the field specs requested by the material model and
     *  the contact model, as well as those provided in the fieldSpecs input argument.
     */
    void initializeDataManager(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs);

    /*! \brief Import data from the given vector to the underlying Epetra_Vector associated with the given field spec.
     *
     *  The intended use case is to import from a non-overlapped vector (i.e., no ghosts) to an overlap vector in
     *  the Block's DataManager.  Note that if the Block does not have space allocated for the given spec and step,
     *  then the function is a no-op.
     */
    void import(const Epetra_Vector& source, Field_NS::FieldSpec spec, Field_ENUM::Step step, Epetra_CombineMode combineMode);

  protected:
    
    std::string blockName;

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

    //! The DataManager
    Teuchos::RCP<PeridigmNS::DataManager> dataManager;

    //! The material model
    Teuchos::RCP<const PeridigmNS::Material> materialModel;

    //! The contact model
    Teuchos::RCP<const PeridigmNS::ContactModel> contactModel;
  };
}

#endif // PERIDIGM_BLOCK_HPP
