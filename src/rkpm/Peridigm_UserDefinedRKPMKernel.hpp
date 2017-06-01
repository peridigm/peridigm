/*! \file Peridigm_UserDefinedRKPMKernel.hpp */

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

#ifndef PERIDIGM_USERDEFINEDRKPMKERNEL_HPP
#define PERIDIGM_USERDEFINEDRKPMKERNEL_HPP

#include "Peridigm_RKPMKernel.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_assert.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include <Trilinos_version.h>
#if TRILINOS_MAJOR_MINOR_VERSION >= 11100
#include "RTC_FunctionRTC.hh"
#else
#include "FunctionRTC.hh"
#endif

namespace PeridigmNS {

  class Peridigm;

  //! Base class defining the Peridigm RKPM Kernel interface.
  class UserDefinedRKPMKernel : public RKPMKernel{

  public:
	
    //! Standard constructor.
    UserDefinedRKPMKernel(const Teuchos::ParameterList& params);

    //! Destructor.
    virtual ~UserDefinedRKPMKernel();

    //! Return name of the model.
    virtual std::string Name() const { return("User Defined"); }

    //! Returns a vector of field IDs corresponding to the variables associated with the model.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Initialize the RKPM Kernel.
    virtual void
    initialize(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               PeridigmNS::DataManager& dataManager) const ;

    //! Evaluate the shape function
    virtual void
    computeRKPMShapeFunction(const double dt,
                             const int numOwnedPoints,
                             const int* ownedIDs,
                             const int* neighborhoodList,
                             PeridigmNS::DataManager& dataManager) const ;

    virtual void
	applyRKPMShapeFunction(const double dt,
                           const int numOwnedPoints,
                           const int* ownedIDs,
                           const int* neighborhoodList,
                           PeridigmNS::DataManager& dataManager) const;

    virtual void
	updateRKPMShapeFunction(const double dt,
                            const int numOwnedPoints,
                            const int* ownedIDs,
                            const int* neighborhoodList,
                            PeridigmNS::DataManager& dataManager) const;

    //! evaluate Parser
    void evaluateParserKernel();

  protected:
    //! Function for evaluating user-defined RKPM Kernel
    static double userdefinedkernel(const double zeta);

    double m_RKPMSupport;
    int m_RKPMBasisOrder;

    //! string defined function;
    std::string functionkernel;

    //! Run-time compiler, used as function parser
    static PG_RuntimeCompiler::Function rtcFunction;

    //! make sure rtcFunction is defined only once
    static int CountKernelDefined ;
    
    // field ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_modelCoordinatesFieldId;
    int m_coordinatesFieldId;
    int m_damageFieldId;
    int m_bondDamageFieldId;
    int m_OwnRKPMShapeFunctionFieldId;
    int m_BondRKPMShapeFunctionFieldId;

  };

}

#endif // PERIDIGM_USERDEFINEDRKPMKERNEL_HPP
