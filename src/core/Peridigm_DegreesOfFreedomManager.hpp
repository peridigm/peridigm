/*! \file Peridigm_DegreesOfFreedom.hpp */

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

#ifndef PERIDIGM_DEGREESOFFREEDOM_HPP
#define PERIDIGM_DEGREESOFFREEDOM_HPP

#include <Teuchos_ParameterList.hpp>

namespace PeridigmNS {

//! Singleton class for tracking the degrees of freedom for multiphysics, etc.
class DegreesOfFreedomManager {

public:

  //! Destructor.
  ~DegreesOfFreedomManager(){}

  //! Singleton.
  static DegreesOfFreedomManager& self();

  void initialize(Teuchos::ParameterList& solverParams);

  bool displacementTreatedAsUnknown() { return displacementTreatedAsUnknown_; }

  bool temperatureTreatedAsUnknown() { return temperatureTreatedAsUnknown_; }

  bool concentrationTreatedAsUnknown() { return concentrationTreatedAsUnknown_; }

  bool pressureTreatedAsUnknown() { return pressureTreatedAsUnknown_; }

  int numberOfDisplacementDegreesOfFreedom() { return number_of_displacement_degrees_of_freedom_; }

  int numberOfTemperatureDegreesOfFreedom() { return number_of_temperature_degrees_of_freedom_; }

  int numberOfConcentrationDegreesOfFreedom() { return number_of_concentration_degrees_of_freedom_; }

  int numberOfPressureDegreesOfFreedom() { return number_of_pressure_degrees_of_freedom_; }

  int totalNumberOfDegreesOfFreedom() { return total_number_of_degrees_of_freedom_; }

  int displacementDofOffset() { return displacement_dof_offset_; }

  int concentrationDofOffset() { return concentration_dof_offset_; }

  int temperatureDofOffset() { return temperature_dof_offset_; }

  int pressureDofOffset() { return pressure_dof_offset_; }

  void print();

private:

  //! Private constructor
  DegreesOfFreedomManager()
    : displacementTreatedAsUnknown_(true), temperatureTreatedAsUnknown_(false),
      concentrationTreatedAsUnknown_(false), pressureTreatedAsUnknown_(false),
      number_of_displacement_degrees_of_freedom_(0), number_of_temperature_degrees_of_freedom_(0),
      number_of_concentration_degrees_of_freedom_(0), number_of_pressure_degrees_of_freedom_(0),
      total_number_of_degrees_of_freedom_(0), displacement_dof_offset_(0),
      concentration_dof_offset_(0), temperature_dof_offset_(0), pressure_dof_offset_(0),
      verbose_(false)
  {}

  //! Private and unimplemented to prevent use

  DegreesOfFreedomManager(const DegreesOfFreedomManager&);

  DegreesOfFreedomManager& operator=(const DegreesOfFreedomManager&);

protected:

  bool displacementTreatedAsUnknown_;
  bool temperatureTreatedAsUnknown_;
  bool concentrationTreatedAsUnknown_;
  bool pressureTreatedAsUnknown_;
  int number_of_displacement_degrees_of_freedom_;
  int number_of_temperature_degrees_of_freedom_;
  int number_of_concentration_degrees_of_freedom_;
  int number_of_pressure_degrees_of_freedom_;
  int total_number_of_degrees_of_freedom_;
  int displacement_dof_offset_;
  int concentration_dof_offset_;
  int temperature_dof_offset_;
  int pressure_dof_offset_;
  bool verbose_;
};

}

#endif // PERIDIGM_DEGREESOFFREEDOM_HPP
