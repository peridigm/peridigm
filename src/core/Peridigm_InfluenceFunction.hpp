/*! \file Peridigm_InfluenceFunction.hpp */

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

#ifndef PERIDIGM_INFLUENCEFUNCTION_HPP
#define PERIDIGM_INFLUENCEFUNCTION_HPP

#include <Teuchos_Assert.hpp>
#include <string>

namespace PeridigmNS{

namespace PeridigmInfluenceFunction{

// Built-in influence functions should be implemented here
// and associated with a string in setInfluenceFunction(), below.

static double one(double zeta, double horizon){ return 1.0; }

}

class InfluenceFunction {

public:

  //! Type definition for the function pointer to an influence function
  typedef double (*functionPointer)(double, double);

  //! Singleton.
  static InfluenceFunction & self() {
    static InfluenceFunction influenceFunction;
    return influenceFunction;
  }

  //! Sets the influence function based on the provided string.
  void setInfluenceFunction(std::string influenceFunctionString) {
    if(influenceFunctionString == "One"){
      m_influenceFunction = &PeridigmInfluenceFunction::one;
    }
    else{
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  InfluenceFunction::setInfluenceFunction(), invalid influence function, must be \"One\"\n");
    }
  }

  //! Returns a function pointer to the influence function.
  functionPointer getInfluenceFunction() {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(m_influenceFunction == NULL,
                                "**** Error:  InfluenceFunction::getInfluenceFunction() called prior to calling InfluenceFunction::setInfluenceFunction().\n")
    return m_influenceFunction;
  }

private:

  //! Constructor, private to prevent use (singleton class).
  InfluenceFunction() : m_influenceFunction(NULL) {
    setInfluenceFunction("One");
  }

  //! Private and unimplemented to prevent use
  InfluenceFunction( const InfluenceFunction & );

  //! Private and unimplemented to prevent use
  InfluenceFunction & operator= ( const InfluenceFunction & );

  //! Function pointer to the influence function with the signature:  double function(double zeta, double horizon).
  functionPointer m_influenceFunction;
};

}

#endif // PERIDIGM_INFLUENCEFUNCTION_HPP
