/*! \file Peridigm_DegreesOfFreedom.cpp */

#include "Peridigm_DegreesOfFreedomManager.hpp"

PeridigmNS::DegreesOfFreedomManager& PeridigmNS::DegreesOfFreedomManager::self() {
  static DegreesOfFreedomManager dofManager;
  return dofManager;
}

void PeridigmNS::DegreesOfFreedomManager::initialize(Teuchos::ParameterList& solverParams) {

  displacementTreatedAsUnknown_ = solverParams.get<bool>("Solve For Displacement", true);
  temperatureTreatedAsUnknown_ = solverParams.get<bool>("Solve For Temperature", false);
  pressureTreatedAsUnknown_ = solverParams.get<bool>("Solve For Pressure", false);

  total_number_of_degrees_of_freedom_ = 0;
  if (displacementTreatedAsUnknown_) {
    displacement_dof_offset_ = total_number_of_degrees_of_freedom_;
    total_number_of_degrees_of_freedom_ += 3;
  }
  if (pressureTreatedAsUnknown_) {
    pressure_dof_offset_ = total_number_of_degrees_of_freedom_;
    total_number_of_degrees_of_freedom_ += 1;
  }
  if (temperatureTreatedAsUnknown_) {
    temperature_dof_offset_ = total_number_of_degrees_of_freedom_;
    total_number_of_degrees_of_freedom_ += 1;
  }
}
