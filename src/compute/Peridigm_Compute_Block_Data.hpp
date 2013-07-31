/*! \file Peridigm_Compute_Block_Data.hpp */

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

#ifdef COMPUTE_CLASS

ComputeClass(Block_Data,Compute_Block_Data)

#else

#ifndef PERIDIGM_COMPUTE_BLOCK_DATA_HPP
#define PERIDIGM_COMPUTE_BLOCK_DATA_HPP

#include "Peridigm_Compute.hpp"

namespace PeridigmNS {

  //! Class for tracking the data at a specified point in space
  class Compute_Block_Data : public PeridigmNS::Compute {

  public:
	
    //! Standard constructor.
    Compute_Block_Data( Teuchos::RCP<const Teuchos::ParameterList> params,
                             Teuchos::RCP<const Epetra_Comm> epetraComm_,
                             Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_);

    //! Destructor.
    ~Compute_Block_Data() {}

    //! Returns a vector of field IDs corresponding to the variables associated with the compute class.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Initialize the compute class
    virtual void initialize( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks );

    //! Perform computation
    virtual int compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const;

  private:

    //! Name of variable to be tracked
    std::string m_variable;
    int m_variableLength;
    bool m_variableIsStated;
    std::string m_blockName;
    int m_blockId;

    enum CALCULATION_TYPE {
      UNDEFINED_CALCULATION = 0,
      MINIMUM,
      MAXIMUM,
      SUM
    } m_calculationType ;

    //! Label for output variable
    std::string m_outputLabel;

    //! Field ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_variableFieldId;
    int m_outputFieldId;
  };
}

#endif // PERIDIGM_COMPUTE_BLOCK_DATA_HPP
#endif // COMPUTE_CLASS
