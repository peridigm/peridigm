/*! \file Peridigm_Factory.cpp */

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

#include "Peridigm_Factory.hpp"
#include <Teuchos_XMLParameterListHelpers.hpp>

#ifdef USE_YAML
  #include <Teuchos_YamlParser_decl.hpp>
  #include <Teuchos_YamlParameterListCoreHelpers.hpp>
#endif

Teuchos::RCP<PeridigmNS::Peridigm> PeridigmNS::PeridigmFactory::create(const std::string inputFile,
                                                                       const MPI_Comm& comm,
                                                                       Teuchos::RCP<Discretization> inputPeridigmDiscretization)
{
  // Input files are read into a ParameterList
  Teuchos::RCP<Teuchos::ParameterList> peridigmParams = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::Ptr<Teuchos::ParameterList> peridigmParamsPtr(peridigmParams.get());
  setPeridigmParamDefaults(peridigmParams.ptr());

  std::string file_extension;
  size_t pos = inputFile.rfind(".");
  if(pos != std::string::npos) {
    file_extension = inputFile.substr(pos);
  }

  if(file_extension == ".xml") {
    Teuchos::updateParametersFromXmlFile(inputFile, peridigmParamsPtr);
  }
  else if(file_extension == ".yaml") {
#ifdef USE_YAML
    Teuchos::updateParametersFromYamlFile(inputFile, peridigmParamsPtr);
#else
    std::string msg = "**** Error, YAML reader not available.  Trilinos must be compiled with YAML support to enable this feature.\n";
                msg += "**** You must obtain the yaml-cpp library, build it, and re-compile Trilinos with the the following flags:\n";
                msg += "****   TPL_ENABLE_yaml-cpp:BOOL=ON\n";
                msg += "****   yaml-cpp_INCLUDE_DIRS:PATH=/path/to/yaml/include\n";
                msg += "****   yaml-cpp_LIBRARY_DIRS:PATH=/path/to/yaml/lib\n";
    throw std::runtime_error(msg);
#endif
  }
  else if(file_extension == ".peridigm") {
    std::string msg = "**** Error, the .peridigm file format is no longer supported.\n**** Input files in the .peridigm format can be converted to .yaml files using the peridigm_to_yaml.py script in the scripts directory.\n";
    throw std::runtime_error(msg);
  }
  else {
    std::string msg = "**** Error, unrecognized input file extension.  Supported formats are .xml and .yaml\n";
    throw std::runtime_error(msg);
  }

  // Create the Peridigm object using the ParameterList
  return Teuchos::rcp(new PeridigmNS::Peridigm(comm, peridigmParams, inputPeridigmDiscretization));
}

Teuchos::RCP<PeridigmNS::Peridigm> PeridigmNS::PeridigmFactory::create(const std::string inputFile,
                                                                       const MPI_Comm& comm)
{
  Teuchos::RCP<Discretization> nullDiscretization;
  return create(inputFile, comm, nullDiscretization);
}

void PeridigmNS::PeridigmFactory::setPeridigmParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> peridigmParams_)
{
  peridigmParams_->set("Verbose", false);
}
