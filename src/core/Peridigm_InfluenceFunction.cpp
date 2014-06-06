/*! \file Peridigm_InfluenceFunction.cpp */

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

#include "Peridigm_InfluenceFunction.hpp"

using namespace std;

// double PeridigmNS::InfluenceFunction::muParserZeta = 0.0;
// double PeridigmNS::InfluenceFunction::muParserHorizon = 0.0;
// mu::Parser PeridigmNS::InfluenceFunction::muParser;

PG_RuntimeCompiler::Function PeridigmNS::InfluenceFunction::rtcFunction(3, "rtcInfluenceFunctionUserDefinedFunction");

PeridigmNS::InfluenceFunction& PeridigmNS::InfluenceFunction::self() {
  static InfluenceFunction influenceFunction;
  return influenceFunction;
}

PeridigmNS::InfluenceFunction::InfluenceFunction() : m_influenceFunction(NULL) {

  // Set the influence function to One by default
  setInfluenceFunction("One");

  // Set up the muParser
  // try {
  //   muParser.DefineVar("zeta", &muParserZeta);
  //   muParser.DefineVar("horizon", &muParserHorizon);
  // } 
  // catch (mu::Parser::exception_type &e)
  //   TEUCHOS_TEST_FOR_EXCEPT_MSG(1, e.GetMsg());

  // set up RTCompiler
  rtcFunction.addVar("double", "zeta");
  rtcFunction.addVar("double", "horizon");
  rtcFunction.addVar("double", "value");
}

double PeridigmNS::InfluenceFunction::userDefinedInfluenceFunction(double zeta, double horizon){
  // muParserZeta = zeta;
  // muParserHorizon = horizon;
  // double value;
  // try {
  //   value = muParser.Eval();
  // }
  // catch (mu::Parser::exception_type &e)
  //   TEUCHOS_TEST_FOR_EXCEPT_MSG(true, e.GetMsg());

  double value(0.0);
  bool success = rtcFunction.varValueFill(0, zeta);
  if(success)
    success = rtcFunction.varValueFill(1, horizon);
  if(success)
    success = rtcFunction.execute();
  if(success)
    value = rtcFunction.getValueOfVar("value");
  return value;
}
