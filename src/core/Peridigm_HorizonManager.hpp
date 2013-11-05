/*! \file Peridigm_Horizon.hpp */

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

#ifndef PERIDIGM_HORIZONMANAGER_HPP
#define PERIDIGM_HORIZONMANAGER_HPP

#include <Teuchos_ParameterList.hpp>
#include "muParser/muParser.h"
#include "muParser/muParserPeridigmFunctions.h"
#include <string>
#include <map>

namespace PeridigmNS{

class HorizonManager {

public:

  //! Constructor
  HorizonManager() : muParserX(0.0), muParserY(0.0), muParserZ(0.0)
  {
    // Set up muParser
    try {
      muParser.DefineVar("x", &muParserX);
      muParser.DefineVar("y", &muParserY);
      muParser.DefineVar("z", &muParserZ);
    }
    catch (mu::Parser::exception_type &e)
      TEUCHOS_TEST_FOR_EXCEPT_MSG(1, e.GetMsg());
  }

  //! Singleton.
  static HorizonManager & self();

  //! Parse the block parameter list and record the horizon information for each block.
  void loadHorizonInformationFromBlockParameters(Teuchos::ParameterList& blockParams);

  //! Returns true if the block has a single, constant value for the horizon.
  bool blockHasConstantHorizon(std::string blockName);

  //! Returns the single, constant horizon value for the block if there is one; throws an exception if the horizon for the given block is not constant.
  double getBlockConstantHorizonValue(std::string blockName);

  //! Evaluates the horizon for a given block at the given coordinates (x, y, z).
  double evaluateHorizon(std::string blockName, double x, double y, double z);

  //! Throws a warning if it seems like the horizon is too big
  // void checkHorizon(Teuchos::RCP<Discretization> peridigmDisc, std::map<std::string, double> & blockHorizonValues);

protected:

  //! Function parser
  mu::Parser muParser;

  //! @name Variables for function parser.
  //@{
  double muParserX;
  double muParserY;
  double muParserZ;
  double muParserT;
  //@}

  //! Container for strings defining horizon for each block.
  std::map<std::string, std::string> horizonStrings;

  //! Record of which blocks have constant horizons.
  std::map<std::string, bool> horizonIsConstant;
};

}

#endif // PERIDIGM_HORIZONMANAGER_HPP
