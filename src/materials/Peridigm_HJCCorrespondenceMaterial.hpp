// SPDX-License-Identifier: BSD-3-Clause
#ifndef PERIDIGM_HJCCORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_HJCCORRESPONDENCEMATERIAL_HPP

#include "Peridigm_CorrespondenceMaterial.hpp"
#include "hjc.hpp"

namespace PeridigmNS {

class HJCCorrespondenceMaterial : public CorrespondenceMaterial {
public:
  explicit HJCCorrespondenceMaterial(const Teuchos::ParameterList& params);
  std::string Name() const override { return "HJC Correspondence"; }
  void initialize(double dt, int numOwnedPoints, const int* ownedIDs,
                  const int* neighborhoodList, DataManager& dataManager) override;
  void computeCauchyStress(double dt, int numOwnedPoints,
                          DataManager& dataManager) const override;
private:
  // These members are used by the out-of-line methods in the .cpp file.
  // cppcheck-suppress unusedStructMember
  hjc::Parameters m_parameters;
  // log volume, plastic volume, maximum compression, density measure, D, EPSP.
  // cppcheck-suppress unusedStructMember
  int m_history[6];
  // cppcheck-suppress unusedStructMember
  int m_pressure, m_equivalentStress, m_acousticModulus;
};

}
#endif
