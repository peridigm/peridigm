// SPDX-License-Identifier: BSD-3-Clause
#include "Peridigm_HJCCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "HjcTests.h"
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_UnitTestRepository.hpp>
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>

using namespace PeridigmNS;

static Teuchos::ParameterList parameters() {
  Teuchos::ParameterList p;
  p.set("Density",2700.); p.set("Shear Modulus",2.417e10);
  p.set("Hourglass Coefficient",.05);
  p.set("A",.29);p.set("B",2.06);p.set("C",.0013);p.set("N",.866);
  p.set("Compressive Strength",1.19e8);p.set("Tensile Strength",8.2e6);
  p.set("Reference Strain Rate",1.);p.set("EFMIN",.01);p.set("SFMAX",5.);
  p.set("Crushing Pressure",4.e7);p.set("Crushing Volumetric Strain",.00124);
  p.set("Locking Pressure",1.2e9);p.set("Locking Plastic Volumetric Strain",.011);
  p.set("D1",4.);p.set("D2",1.);
  p.set("K1",1.287e10);p.set("K2",1.631e10);p.set("K3",6.495e10);
  return p;
}

TEUCHOS_UNIT_TEST(HJC, analyticBranches) {
  TEST_NOTHROW(hjc_tests::runAnalytic());
}

TEUCHOS_UNIT_TEST(HJC, parameterValidation) {
  Teuchos::ParameterList p=parameters();
  TEST_NOTHROW(HJCCorrespondenceMaterial{p});
  p.set("Crushing Volumetric Strain",0.);
  TEST_THROW(HJCCorrespondenceMaterial{p},std::exception);
  p=parameters();p.set("Bulk Modulus",1.);
  TEST_THROW(HJCCorrespondenceMaterial{p},std::exception);
  p=parameters();p.set("Thermal Expansion Coefficient",1.e-5);
  TEST_THROW(HJCCorrespondenceMaterial{p},std::exception);
  p=parameters();p.set("EFMIN",0.);
  TEST_THROW(HJCCorrespondenceMaterial{p},std::exception);
  p=parameters();p.set("Crushing Volumetric Strain",.1);
  TEST_THROW(HJCCorrespondenceMaterial{p},std::exception);
}

TEUCHOS_UNIT_TEST(HJC, historyAndTensorMapping) {
  HJCCorrespondenceMaterial material(parameters());
  Epetra_SerialComm comm;
  Epetra_Map nodes(8,0,comm),vectors(24,0,comm),bonds(56,0,comm);
  DataManager dm;
  dm.setMaps(Teuchos::rcp(&nodes,false),Teuchos::rcp(&nodes,false),
    Teuchos::rcp(&vectors,false),Teuchos::rcp(&vectors,false),Teuchos::rcp(&bonds,false));
  dm.allocateData(material.FieldIds());
  FieldManager& fm=FieldManager::self();
  auto data=[&](const char* name,PeridigmField::Step step)->Epetra_Vector& {
    return *dm.getData(fm.getFieldId(name),step);
  };
  auto& x=data("Model_Coordinates",PeridigmField::STEP_NONE);
  auto& y=data("Coordinates",PeridigmField::STEP_NP1);
  for(int i=0;i<8;++i) for(int axis=0;axis<3;++axis)
    y[3*i+axis]=x[3*i+axis]=(i>>axis)&1;
  data("Volume",PeridigmField::STEP_NONE).PutScalar(1.);
  data("Horizon",PeridigmField::STEP_NONE).PutScalar(2.);
  int owned[8],neighbors[64];
  for(int i=0;i<8;++i) {
    owned[i]=i;neighbors[8*i]=7;int k=8*i+1;
    for(int j=0;j<8;++j) if(j!=i) neighbors[k++]=j;
  }
  material.initialize(1.e-6,8,owned,neighbors,dm);
  const double initialDamage=0.;
  TEST_FLOATING_EQUALITY(data("HJC_Damage",PeridigmField::STEP_N)[0],initialDamage,1.e-13);
  auto& rate=data("Unrotated_Rate_Of_Deformation",PeridigmField::STEP_NONE);
  for(int i=0;i<8;++i) {
    rate[9*i+1]=rate[9*i+3]=.5;
    rate[9*i+2]=rate[9*i+6]=-.25;
    rate[9*i+5]=rate[9*i+7]=.75;
  }
  material.computeCauchyStress(1.e-6,8,dm);
  auto& stress=data("Unrotated_Cauchy_Stress",PeridigmField::STEP_NP1);
  TEST_FLOATING_EQUALITY(stress[1],2.417e4,1.e-12);
  TEST_FLOATING_EQUALITY(stress[2],-1.2085e4,1.e-12);
  TEST_FLOATING_EQUALITY(stress[5],3.6255e4,1.e-12);
  TEST_EQUALITY(stress[1],stress[3]);TEST_EQUALITY(stress[2],stress[6]);
  material.computeCauchyStress(1.e-6,8,dm);
  TEST_FLOATING_EQUALITY(stress[1],2.417e4,1.e-12); // trial evaluation is repeatable
  TEST_EQUALITY(data("Equivalent_Plastic_Strain",PeridigmField::STEP_N)[0],0.);
  dm.updateState();
  rate.PutScalar(0.);
  material.computeCauchyStress(0.,8,dm);
  TEST_EQUALITY(data("Unrotated_Cauchy_Stress",PeridigmField::STEP_NP1)[1],2.417e4);
  TEST_EQUALITY(data("HJC_Damage",PeridigmField::STEP_NP1)[0],initialDamage);
  TEST_THROW(material.computeCauchyStress(-1.,8,dm),std::exception);
}

int main(int argc,char** argv) {
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc,argv);
}
