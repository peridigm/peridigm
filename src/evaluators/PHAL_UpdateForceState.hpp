// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
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
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PHAL_UPDATEFORCESTATE_HPP
#define PHAL_UPDATEFORCESTATE_HPP

#include <vector>

#include <Phalanx_ConfigDefs.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>
#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_MDField.hpp>
#include "PHAL_ParameterEntry.hpp"
#include "PHAL_ParameterGet.hpp"

template<typename EvalT, typename Traits>
class UpdateForceState : public PHX::EvaluatorWithBaseImpl<Traits>,
						 public PHX::EvaluatorDerived<EvalT, Traits>,
						 public PHAL::ParameterGet<EvalT>
{
  typedef typename EvalT::ScalarT ScalarT;

public:

  UpdateForceState(const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename Traits::SetupData d,
							 PHX::FieldManager<Traits>& vm);

  void evaluateFields(typename Traits::EvalData d);

  // getValue method is pure virtual in PHAL::ParameterGet
  virtual ScalarT & getValue(const std::string &n) {
	return parameter_accessible_via_getValue;
  }

private:

  ScalarT parameter_accessible_via_getValue;
  bool m_verbose;
 
  Teuchos::RCP<PHX::FieldTag> update_force_state_field_tag;

  std::size_t m_num_pt;

  void setup_vectors(const Teuchos::ParameterList& p);
};

#include "PHAL_UpdateForceState_Def.hpp"

#endif // PHAL_UPDATEFORCESTATE_HPP
