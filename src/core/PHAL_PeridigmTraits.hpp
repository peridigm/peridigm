/*! \file PHAL_PeridigmTraits.hpp */

// @HEADER
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
#ifndef PHAL_PERIDIGMTRAITS_HPP
#define PHAL_PERIDIGMTRAITS_HPP

// mpl (Meta Programming Library) templates
#include <Sacado.hpp>
#include <Sacado_mpl_vector.hpp>
#include <Sacado_mpl_find.hpp>
#include <boost/mpl/map.hpp>
//#include <boost/mpl/find.hpp>

// traits Base Class
#include <Phalanx_Traits_Base.hpp>

// Include User Data Types
//#include <Phalanx_ConfigDefs.hpp>
#include <Phalanx_Allocator_Contiguous.hpp>
#include <Phalanx_TypeStrings.hpp>

//#include "PHAL_Dimension.hpp"
#include "PHAL_Workset.hpp"

// Include all of our AD types
//#include <Sacado_Fad_DFad.hpp>

// Typedef AD types to standard names
typedef double RealType;
typedef Sacado::Fad::DFad<double> FadType;

// Include ScalarParameterLibrary to specialize its traits
#include <Sacado_ScalarParameterLibrary.hpp>

namespace PHAL {

  /*! The Traits object is a struct that defines the evaluation types, 
   *  the data types for each evaluation type, the allocator type, and the 
   *  user defined types that will be passed through evaluator calls.
   */
  struct PeridigmTraits : public PHX::TraitsBase {
    
    // ******************************************************************
    // *** Evaluation Types
    // ******************************************************************

    // Create the scalar types for each evaluation type

    struct Residual { typedef RealType ScalarT; };
    struct Jacobian { typedef FadType  ScalarT; };
    struct Tangent  { typedef FadType  ScalarT; };

   /*! EvalTypes - an mpl::vector of user defined evaluation types. Each 
    *  evaluation type must have a typedef member called ScalarT that provides 
    *  the default scalar type. This is used to automate the building of evaluators
    *  for each evaluation type using the EvaluatorFactory.
    */
    typedef Sacado::mpl::vector<Residual> EvalTypes;

    // In general, we'd want to define types for residual, jacobian, and tagent
//    typedef Sacado::mpl::vector<Residual, Jacobian, Tangent> EvalTypes;

    // ******************************************************************
    // *** Eval to Data Map
    // ******************************************************************
    
    // Create the data types for each evaluation type
    
    // Residual (default scalar type is RealType)
    typedef Sacado::mpl::vector<RealType> ResidualDataTypes;
  
    // Jacobian (default scalar type is Fad<double, double>)
    typedef Sacado::mpl::vector<FadType> JacobianDataTypes;

    // Tangent (default scalar type is Fad<double, double>)
    typedef Sacado::mpl::vector<FadType> TangentDataTypes;

    // Maps the key EvalType a vector of DataTypes

   /*! EvalToDataMap - an mpl::map. The key is an evaluation type and the value 
    *  is an mpl::vector of valid data types for that particular evaluation type.
    */
    typedef boost::mpl::map< boost::mpl::pair<Residual, ResidualDataTypes> >::type EvalToDataMap;

    // In general, we'd want to define types for residual, jacobian, and tagent
//    typedef boost::mpl::map<
//      boost::mpl::pair<Residual, ResidualDataTypes>,
//      boost::mpl::pair<Jacobian, JacobianDataTypes>, 
//      boost::mpl::pair<Tangent, TangentDataTypes>
//    >::type EvalToDataMap;

    // ******************************************************************
    // *** Allocator Type
    // ******************************************************************

    /*! Allocator type - type that defines the allocator class to use to allocate
     *  the memory for data storage.  Phalanx comes with two allocators, but the user
     *  can write their own allocator class if these aren't sufficient. The Phalanx 
     *  allocator classes are:
     *  PHX::NewAllocator: uses the C++ "new" command to allocate each field on the heap 
     *  separately.
     *  PHX::ContiguousAllocator: allocates a single contiguous block of memory on the heap 
     *  for all fields regardless of the type. This allows us to fit a subset of elements 
     *  into cache to speed up the evaluation.
     */
    typedef PHX::ContiguousAllocator<double> Allocator;

    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************

   /*! Users can pass their own data to the evaluateFields(), preEvaluate()
    *  and postEvaluate() methods of the PHX::FiledManager class. In this 
    *  example, the user passes in a struct that they have written called 
    *  MyEvalData. This contains information about the cell workset. The user is
    *  not required to write their own object. they could just pass in a null 
    *  pointer if they don't need auxiliary information passed into the routine. 
    *  This is demonstrated in the PreEvalData and PostEvalData. A void* is set
    *  for the data member.
    */

   /*! EvalData - A user defined type to be passed in to the evaluateFields()
    *  call. Allows users to pass in arbitrary data.
    */
    typedef const Workset& EvalData;

   /*! PreEvalData - A user defined type to be passed in to the preEvaluate() 
    *  call. Allows users to pass in arbitrary data.
    */
    typedef void* PreEvalData;

   /*! PostEvalData - A user defined type to be passed in to the postEvaluate() 
    *  call. Allows users to pass in arbitrary data.
    */
    typedef void* PostEvalData;

   /*! SetupData - A user defined type to be passed in to the postRegistrationSetup() 
    *  call. Allows users to pass in arbitrary data.
    */
    typedef void* SetupData;
  };
 
  // ******************************************************************
  // ******************************************************************
  // Debug strings.  Specialize the Evaluation and Data types for the
  // TypeString object in phalanx/src/Phalanx_TypeStrings.hpp.
  // ******************************************************************
  // ******************************************************************

}

namespace PHX {

   /*! Specialize the PHX::TypeString Object
    *  For debugging information, Phalanx makes a forward declaration of the
    *  PHX::TypeString object. This must be specialized for each evaluation type
    *  and each data type so that if there is a run-time error, phalanx can report 
    *  detailed information on the problem. We could have used the typeinfo from the
    *  stl, but the name() method is not demangled on every platform, so it can make 
    *  debugging a challenge. The specialized classes can go into their own file or 
    *  can be added to the traits file above depending on how you use the Traits class.
    *  During linking, if the compiler complains about multiple defninitions of your 
    *  specialized traits classes, separate the traits implementation into their own file.
    */

  // Evaluation Types
  template<> struct TypeString<PHAL::PeridigmTraits::Residual> 
  { static const std::string value; };

  template<> struct TypeString<PHAL::PeridigmTraits::Jacobian> 
  { static const std::string value; };

  template<> struct TypeString<PHAL::PeridigmTraits::Tangent> 
  { static const std::string value; };

  // Data Types
  template<> struct TypeString<double> 
  { static const std::string value; };

  template<> struct TypeString< Sacado::Fad::DFad<double> > 
  { static const std::string value; };

 /*! String values are defined in PHAL_PeridigmTraits.cpp. */

}

// ******************************************************************
// Definition of Sacado::ParameterLibrary traits
// ******************************************************************
struct SPL_Traits {
  template <class T> struct apply {
    typedef typename T::ScalarT type;
  };
};

// Synonym for the ScalarParameterLibrary/Vector on our traits
typedef Sacado::ScalarParameterLibrary<SPL_Traits> ParamLib;
typedef Sacado::ScalarParameterVector<SPL_Traits> ParamVec;

// Turn on/off explicit template instantiation
#define SACADO_ETI

// Define macro for explicit template instantiation
#define INSTANTIATE_TEMPLATE_CLASS_RESIDUAL(name) template class name<PHAL::PeridigmTraits::Residual>;
#define INSTANTIATE_TEMPLATE_CLASS_JACOBIAN(name) template class name<PHAL::PeridigmTraits::Jacobian>;

#define INSTANTIATE_TEMPLATE_CLASS(name) \
  INSTANTIATE_TEMPLATE_CLASS_RESIDUAL(name)	 \
  INSTANTIATE_TEMPLATE_CLASS_JACOBIAN(name)

#endif
