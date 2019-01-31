/*! \file Peridigm_DegreesOfFreedom.cpp */

#include "Peridigm_DegreesOfFreedomManager.hpp"
#include <iostream>

PeridigmNS::DegreesOfFreedomManager& PeridigmNS::DegreesOfFreedomManager::self() {
  static DegreesOfFreedomManager dofManager;
  return dofManager;
}

void PeridigmNS::DegreesOfFreedomManager::initialize(Teuchos::ParameterList& solverParams) {

  displacementTreatedAsUnknown_ = solverParams.get<bool>("Solve For Displacement", true);
  temperatureTreatedAsUnknown_ = solverParams.get<bool>("Solve For Temperature", false);
  pressureTreatedAsUnknown_ = solverParams.get<bool>("Solve For Pressure", false);

  if (displacementTreatedAsUnknown_) {
    displacement_dof_offset_ = total_number_of_degrees_of_freedom_;
    number_of_displacement_degrees_of_freedom_ = 3;
    total_number_of_degrees_of_freedom_ += 3;
  }
  if (temperatureTreatedAsUnknown_) {
    temperature_dof_offset_ = total_number_of_degrees_of_freedom_;
    number_of_temperature_degrees_of_freedom_ = 1;
    total_number_of_degrees_of_freedom_ += 1;
  }
  if (pressureTreatedAsUnknown_) {
    pressure_dof_offset_ = total_number_of_degrees_of_freedom_;
    number_of_pressure_degrees_of_freedom_ = 1;
    total_number_of_degrees_of_freedom_ += 1;
  }

  if (solverParams.isParameter("Solve For Displacement") || solverParams.isParameter("Solve For Temperature") || solverParams.isParameter("Solve For Pressure")) {
    verbose_ = true;
  }
}

void PeridigmNS::DegreesOfFreedomManager::print() {
  if (verbose_) {
    std::cout << "Degrees of freedom:" << std::endl;
    std::cout << "  number of displacement degrees of freedom: " << number_of_displacement_degrees_of_freedom_ << std::endl;
    std::cout << "  number of temperature degrees of freedom:  " << number_of_temperature_degrees_of_freedom_ << std::endl;
    std::cout << "  number of pressure degrees of freedom:     " << number_of_pressure_degrees_of_freedom_ << "\n" << std::endl;
  }
}
