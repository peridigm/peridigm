/*! \file PHAL_ParameterEntry.hpp */

// $Source: /space/CVS/Trilinos/packages/sacado/example/FEApp/FEApp_ConstantNodeBCStrategy.hpp,v $ 

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

#include <Sacado_Traits.hpp>
#include <sstream>
#include "PHAL_ParameterGet.hpp"
#include "PHAL_PeridigmTraits.hpp"

namespace PHAL {
  /*!
   * @brief Parameter class for sensitivity/stability analysis 
   */
  template <typename EvalT>
  class ParameterEntry : 
    public Sacado::ScalarParameterEntry<EvalT,SPL_Traits> {

  //! Scalar type
  typedef typename Sacado::ScalarParameterEntry<EvalT,SPL_Traits>::ScalarT ScalarT;

  public:

    //! Constructor
    ParameterEntry(const std::string &name_, ParameterGet<EvalT>* evaluator_,
                            Teuchos::RCP<ParamLib> paramLib)
      : name(name_), evaluator(evaluator_) {

      if (paramLib != Teuchos::null) {
        if (!paramLib->isParameter(name))
          paramLib->addParameterFamily(name, true, false);
        if (!paramLib->template isParameterForType<EvalT>(name)) {
          paramLib->template addEntry<EvalT>(name, Teuchos::rcp(this,false));
        }
      }
    }

    //! Destructor
    virtual ~ParameterEntry() {}

    //! Set real parameter value
    virtual void setRealValue(double value) { 
      setValue(ScalarT(value)); }

    //! Set parameter this object represents to \em value
    virtual void setValue(const ScalarT& value) { 
      evaluator->getValue(name) =  value; }
    
    //! Get parameter value this object represents
    virtual const ScalarT& getValue() const { 
      return evaluator->getValue(name); 
     }
    
  protected:  
    
    //! Pointer to source function
    const std::string name;
    ParameterGet<EvalT>* evaluator;
  };
}

#endif // PHAL_PARAMETERENTRY_HPP
