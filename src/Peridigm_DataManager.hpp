/*! \file Peridigm_DataManager.hpp */

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

#ifndef PERIDIGM_DATAMANAGER_HPP
#define PERIDIGM_DATAMANAGER_HPP

#include "Field.h"
#include "Peridigm_State.hpp"

namespace PeridigmNS {

/*! \brief A lean, mean, data managing machine.
 *
 * DataManager is a container class for managing the data required for peridynamic force evaluations. 
 * Data is stored in State objects which may be either stateless (STATE_NONE) or stateful (STATE_N, STATE_NP1).
 * Three types of data are supported:  scalar, vector, and bond data.  In the cases of scalar and vector data,
 * the DataManger class has notions of owned points and ghosted points, allowing force calcuations access
 * to off-processor points that are within the neighborhood of one or more owned points.
 * Data is accessed via the getData function, which provides access to the Epetra_Vector corresponding to
 * the given FieldSpec and FieldStep.
 */
class DataManager {

public:
  
  //! Constructor.
  DataManager() : rebalanceCount(0) {}

  //! Copy constructor.
  DataManager(const DataManager& dataManager){}

  //! Destructor.
  ~DataManager(){}

  /*! \brief Sets the maps, must be called prior to allocating data.
   *
   * Sets the following maps:
   * <ul>
   * <li> ownedIDScalarMap:  Identifies the points owned by the calling processor.
   * <li> ownedIDVectorMap:  Identical to ownedIDScalarMap, but with dimension three.
   * <li> scalarMap:         One-dimensional overlap map containing the owned points and ghosts.
   * <li> vectorMap:         Three-dimensional overlap map containing the owned points and ghosts.
   * <li> bondMap:           Map containing one element for each owned point; the length of the element is equal to the number of bonds.
   * </ul>
   */
  void setMaps(Teuchos::RCP<const Epetra_BlockMap> ownedIDScalarMap_,
               Teuchos::RCP<const Epetra_BlockMap> ownedIDVectorMap_,
               Teuchos::RCP<const Epetra_BlockMap> scalarMap_,
               Teuchos::RCP<const Epetra_BlockMap> vectorMap_,
               Teuchos::RCP<const Epetra_BlockMap> bondMap_){
    ownedIDScalarMap = ownedIDScalarMap_;
    ownedIDVectorMap = ownedIDVectorMap_;
    scalarMap = scalarMap_;
    vectorMap = vectorMap_;
    bondMap = bondMap_;
  }

  //! Instantiates State objects corresponding to the given list of FieldSpecs. 
  void allocateData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > specs);

  //! For each point, copies data from the processor that owns the point to processors that have ghosted copies.
  void scatterToGhosts();

  //! Redistributes data based on new partitioning as defined by given maps.
  void rebalance(Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedIDScalarMap,
                 Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedIDVectorMap,
                 Teuchos::RCP<const Epetra_BlockMap> rebalancedScalarMap,
                 Teuchos::RCP<const Epetra_BlockMap> rebalancedVectorMap,
                 Teuchos::RCP<const Epetra_BlockMap> rebalancedBondMap);

  //! Returns the number of times rebalance has been called.
  int getRebalanceCount(){ return rebalanceCount; }

  //! Returns RCP to owned scalar map
  Teuchos::RCP<const Epetra_BlockMap> getOwnedIDScalarMap(){ return ownedIDScalarMap; }

  //! Returns RCP to owned vector map
  Teuchos::RCP<const Epetra_BlockMap> getOwnedIDVectorMap(){ return ownedIDVectorMap; }

  //! Returns RCP to overlap scalar map
  Teuchos::RCP<const Epetra_BlockMap> getOverlapIDScalarMap(){ return scalarMap; }

  //! Provides access to the Epetra_Vector specified by the given FieldSped and FieldStep.
  Teuchos::RCP<Epetra_Vector> getData(Field_NS::FieldSpec fieldSpec,
                                      Field_NS::FieldSpec::FieldStep fieldStep);

  //! Returns the complete list of field specs.
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > getFieldSpecs() { return fieldSpecs; }

  //! Swaps StateN and StateNP1; stateNONE is unaffected.
  void updateState(){
    stateN.swap(stateNP1);
  }

protected:

  //! Number of times rebalance has been called.
  int rebalanceCount;

  //! @name Field specifications
  //@{
  //! Complete list of field specs.
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs;
  //! Field specs for stateless scalar data.
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessScalarFieldSpecs;
  //! Field specs for stateless vector data.
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessVectorFieldSpecs;
  //! Field specs for stateless bond data.
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statelessBondFieldSpecs;
  //! Field specs for stateful scalar data.
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulScalarFieldSpecs;
  //! Field specs for stateful vector data.
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulVectorFieldSpecs;
  //! Field specs for stateful bond data.
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > statefulBondFieldSpecs;
  //@}

  //! @name Maps
  //@{
  //! One-dimensional map for owned points.
  Teuchos::RCP<const Epetra_BlockMap> ownedIDScalarMap;
  //! Three-dimensional map for owned points.
  Teuchos::RCP<const Epetra_BlockMap> ownedIDVectorMap;
  //! One-dimensional overlap map for owned points and ghosts.
  Teuchos::RCP<const Epetra_BlockMap> scalarMap;
  //! Three-dimensional overlap map for owned points and ghosts.
  Teuchos::RCP<const Epetra_BlockMap> vectorMap;
  //! Bond map for owned points; the map contains one element for each owned point, the element length is equal to the number of bonds for that point.
  Teuchos::RCP<const Epetra_BlockMap> bondMap;
  //@}

  //! @name State objects
  //@{
  //! Data storage for state N.
  Teuchos::RCP<State> stateN;
  //! Data storage for state N plus 1.
  Teuchos::RCP<State> stateNP1;
  //! Data storage for state NONE (stateless data).
  Teuchos::RCP<State> stateNONE;
  //@}
};

}

#endif // PERIDIGM_DATAMANAGER_HPP
