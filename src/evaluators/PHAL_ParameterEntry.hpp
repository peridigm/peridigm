// $Source: /space/CVS/Trilinos/packages/sacado/example/FEApp/FEApp_ConstantNodeBCStrategy.hpp,v $ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef PHAL_PARAMETERENTRY_HPP
#define PHAL_PARAMETERENTRY_HPP

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
