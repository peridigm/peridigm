// SPDX-License-Identifier: BSD-3-Clause
#include "Peridigm_HJCCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include <Teuchos_Assert.hpp>
#include <cmath>

namespace {
Teuchos::ParameterList elasticParameters(const Teuchos::ParameterList& input) {
  Teuchos::ParameterList p(input);
  const double pc=p.get<double>("Crushing Pressure");
  const double muc=p.get<double>("Crushing Volumetric Strain");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!(pc>0.0 && muc>0.0),
    "HJC requires positive crushing pressure and volumetric strain.");
  const double bulk=pc/muc;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(p.isParameter("Bulk Modulus") &&
    std::fabs(p.get<double>("Bulk Modulus")-bulk)>1.e-8*bulk,
    "HJC Bulk Modulus must equal Crushing Pressure / Crushing Volumetric Strain.");
  p.set("Bulk Modulus",bulk);
  return p;
}
}

PeridigmNS::HJCCorrespondenceMaterial::HJCCorrespondenceMaterial(
    const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(elasticParameters(params)) {
  hjc::Parameters& p=m_parameters;
  p.g=m_shearModulus;
  p.a=params.get<double>("A"); p.b=params.get<double>("B");
  p.c=params.get<double>("C"); p.n=params.get<double>("N");
  p.fc=params.get<double>("Compressive Strength");
  p.tension=params.get<double>("Tensile Strength");
  p.rate0=params.get<double>("Reference Strain Rate");
  p.efmin=params.get<double>("EFMIN"); p.sfmax=params.get<double>("SFMAX");
  p.pc=params.get<double>("Crushing Pressure");
  p.muc=params.get<double>("Crushing Volumetric Strain");
  p.pl=params.get<double>("Locking Pressure");
  p.mul=params.get<double>("Locking Plastic Volumetric Strain");
  p.d1=params.get<double>("D1"); p.d2=params.get<double>("D2");
  p.k1=params.get<double>("K1"); p.k2=params.get<double>("K2");
  p.k3=params.get<double>("K3");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!hjc::prepare(p) ||
    !std::isfinite(m_density) || m_density<=0.0,"Invalid HJC material parameters.");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Thermal Expansion Coefficient"),
    "HJC Correspondence does not include thermal expansion.");

  FieldManager& fm=FieldManager::self();
  const char* names[]={"HJC_Log_Volume","HJC_Plastic_Volume",
    "HJC_Maximum_Compression","HJC_Density_Measure","HJC_Damage",
    "Equivalent_Plastic_Strain"};
  for(int i=0;i<6;++i) {
    m_history[i]=fm.getFieldId(PeridigmField::ELEMENT,PeridigmField::SCALAR,
                             PeridigmField::TWO_STEP,names[i]);
    m_fieldIds.push_back(m_history[i]);
  }
  m_pressure=fm.getFieldId(PeridigmField::ELEMENT,PeridigmField::SCALAR,
                         PeridigmField::CONSTANT,"HJC_Pressure");
  m_equivalentStress=fm.getFieldId(PeridigmField::ELEMENT,PeridigmField::SCALAR,
                                 PeridigmField::CONSTANT,"Von_Mises_Stress");
  m_acousticModulus=fm.getFieldId(PeridigmField::ELEMENT,PeridigmField::SCALAR,
                                PeridigmField::CONSTANT,"HJC_Acoustic_Modulus");
  m_fieldIds.push_back(m_pressure);
  m_fieldIds.push_back(m_equivalentStress);
  m_fieldIds.push_back(m_acousticModulus);
}

void PeridigmNS::HJCCorrespondenceMaterial::initialize(double dt,int n,
    const int* ownedIDs,const int* neighborhoodList,DataManager& dm) {
  CorrespondenceMaterial::initialize(dt,n,ownedIDs,neighborhoodList,dm);
  const hjc::State initial=hjc::initialize(m_parameters);
  for(int i=0;i<6;++i) {
    const double value=i==4 ? initial.damage : 0.0;
    dm.getData(m_history[i],PeridigmField::STEP_N)->PutScalar(value);
    dm.getData(m_history[i],PeridigmField::STEP_NP1)->PutScalar(value);
  }
  dm.getData(m_pressure,PeridigmField::STEP_NONE)->PutScalar(0.0);
  dm.getData(m_equivalentStress,PeridigmField::STEP_NONE)->PutScalar(0.0);
  dm.getData(m_acousticModulus,PeridigmField::STEP_NONE)->PutScalar(
    m_bulkModulus+4.0*m_shearModulus/3.0);
}

void PeridigmNS::HJCCorrespondenceMaterial::computeCauchyStress(double dt,int n,
    DataManager& dm) const {
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(dt) || dt<0.0,
                            "HJC requires a finite nonnegative timestep.");
  double *oldStress,*stress,*rate,*old[6],*next[6],*pressure,*q,*acoustic;
  dm.getData(m_unrotatedCauchyStressFieldId,PeridigmField::STEP_N)->ExtractView(&oldStress);
  dm.getData(m_unrotatedCauchyStressFieldId,PeridigmField::STEP_NP1)->ExtractView(&stress);
  dm.getData(m_unrotatedRateOfDeformationFieldId,PeridigmField::STEP_NONE)->ExtractView(&rate);
  for(int i=0;i<6;++i) {
    dm.getData(m_history[i],PeridigmField::STEP_N)->ExtractView(&old[i]);
    dm.getData(m_history[i],PeridigmField::STEP_NP1)->ExtractView(&next[i]);
  }
  dm.getData(m_pressure,PeridigmField::STEP_NONE)->ExtractView(&pressure);
  dm.getData(m_equivalentStress,PeridigmField::STEP_NONE)->ExtractView(&q);
  dm.getData(m_acousticModulus,PeridigmField::STEP_NONE)->ExtractView(&acoustic);
  const int diagonal[]={0,4,8};
  const int upper[]={1,5,2},lower[]={3,7,6};
  for(int point=0;point<n;++point) {
    hjc::State s={}; double deps[6];
    for(int i=0;i<3;++i) {
      s.stress[i]=oldStress[9*point+diagonal[i]];
      s.stress[i+3]=0.5*(oldStress[9*point+upper[i]]+oldStress[9*point+lower[i]]);
      deps[i]=dt*rate[9*point+diagonal[i]];
      deps[i+3]=dt*(rate[9*point+upper[i]]+rate[9*point+lower[i]]);
    }
    s.volume_strain=old[0][point]; s.plastic_volume=old[1][point];
    s.max_compression=old[2][point]; s.density_measure=old[3][point];
    s.damage=old[4][point]; s.plastic_strain=old[5][point];
    hjc::Result result={};
    if(dt>0.0) result=hjc::update(m_parameters,deps,dt,s);
    else {
      // Initialization/force queries with zero elapsed time must not return-map
      // the old stress at a different strain rate or evolve irreversible state.
      result.pressure=-(s.stress[0]+s.stress[1]+s.stress[2])/3.0;
      double norm=0.0;
      for(int i=0;i<3;++i) norm+=1.5*std::pow(s.stress[i]+result.pressure,2)+3*s.stress[i+3]*s.stress[i+3];
      result.equivalent_stress=std::sqrt(norm);
      result.acoustic_modulus=acoustic[point];
    }
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(result.acoustic_modulus) ||
      result.acoustic_modulus<=0.0 || !std::isfinite(s.volume_strain) ||
      !std::isfinite(s.plastic_strain) || !std::isfinite(s.plastic_volume),
      "Nonfinite HJC history or acoustic modulus; reduce the timestep and check deformation.");
    for(int i=0;i<3;++i) {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(s.stress[i]) ||
        !std::isfinite(s.stress[i+3]),"Nonfinite HJC stress; reduce the timestep and check deformation.");
      stress[9*point+diagonal[i]]=s.stress[i];
      stress[9*point+upper[i]]=stress[9*point+lower[i]]=s.stress[i+3];
    }
    next[0][point]=s.volume_strain; next[1][point]=s.plastic_volume;
    next[2][point]=s.max_compression; next[3][point]=s.density_measure;
    next[4][point]=s.damage; next[5][point]=s.plastic_strain;
    pressure[point]=result.pressure; q[point]=result.equivalent_stress;
    acoustic[point]=result.acoustic_modulus;
  }
}
