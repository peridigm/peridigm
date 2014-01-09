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

#include "Peridigm_BlockBase.hpp"
#include "Peridigm_Material.hpp"
#include "Peridigm_DamageModel.hpp"

namespace PeridigmNS {

  class Block : public BlockBase {
  public:

    //! Constructor
    Block() : BlockBase() {}

    //! Constructor
    Block(std::string blockName_, int blockID_, Teuchos::ParameterList& blockParams_)
      : BlockBase(blockName_, blockID_, blockParams_) {}

    //! Destructor
    ~Block(){}

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

    //! Get the material model
    Teuchos::RCP<const PeridigmNS::Material> getMaterialModel(){
      return materialModel;
    }

    //! Set the material model
    void setMaterialModel(Teuchos::RCP<const PeridigmNS::Material> materialModel_){
      materialModel = materialModel_;
    }

    //! Get the damage model
    Teuchos::RCP<const PeridigmNS::DamageModel> getDamageModel(){
      return damageModel;
    }

    //! Set the damage model
    void setDamageModel(Teuchos::RCP<const PeridigmNS::DamageModel> damageModel_){
      damageModel = damageModel_;
    }

    //! Get the material name
    std::string getMaterialName(){
      return blockParams.get<std::string>("Material");
    }

    //! Get the damage model name
    std::string getDamageModelName(){
      return blockParams.get<std::string>("Damage Model", "None");
    }

    //! Initialize the material model
    void initializeMaterialModel(double timeStep = 1.0);

    //! Initialize the damage model
    void initializeDamageModel(double timeStep = 1.0);

  protected:

    //! The material model
    Teuchos::RCP<const PeridigmNS::Material> materialModel;

    //! The damage model
    Teuchos::RCP<const PeridigmNS::DamageModel> damageModel;
  };
}

#endif // PERIDIGM_BLOCK_HPP
