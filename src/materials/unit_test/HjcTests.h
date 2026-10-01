// SPDX-License-Identifier: BSD-3-Clause
#ifndef HJC_TESTS_H
#define HJC_TESTS_H
#include "hjc.hpp"
#include <limits>
#include <stdexcept>

namespace hjc_tests {
inline hjc::Parameters parameters() {
  hjc::Parameters p={2.417e10,.29,2.06,.0013,.866,1.19e8,8.2e6,1.,.01,5.,
    4.e7,.00124,1.2e9,.011,.04,1.,1.287e10,1.631e10,6.495e10};
  if(!hjc::prepare(p)) throw std::runtime_error("invalid test material");
  return p;
}
inline void check(bool ok,const char* message) {
  if(!ok) throw std::runtime_error(message);
}
inline void close(double a,double b,double atol,const char* message) {
  check(std::isfinite(a) && std::fabs(a-b)<=atol,message);
}

inline void runAnalytic() {
  hjc::Parameters p=parameters();
  check(hjc::valid(p),"valid parameters rejected");
  hjc::Parameters invalid=p; invalid.muc=0;
  check(!hjc::valid(invalid),"zero crushing strain accepted");
  invalid=p; invalid.fc=std::numeric_limits<double>::quiet_NaN();
  check(!hjc::valid(invalid),"NaN accepted");
  hjc::State s=hjc::initialize(p);
  const double hydro[6]={-1.e-6,-1.e-6,-1.e-6,0,0,0};
  hjc::Result r=hjc::update(p,hydro,1.e-6,s);
  close(r.pressure,p.pc/p.muc*std::expm1(3.e-6),1.e-4,"elastic hydrostatic pressure");
  close(s.plastic_volume,0.,0.,"elastic volume must be reversible");
  close(s.plastic_strain,0.,1.e-14,"hydrostatic plastic shear strain");
  s=hjc::initialize(p);
  const double tiny[6]={-1.e-18,-1.e-18,-1.e-18,0,0,0};
  r=hjc::update(p,tiny,1.e-6,s);
  close(r.pressure,3.e-18*p.pc/p.muc,1.e-20,"infinitesimal bulk response");
  for(int axis=3;axis<6;++axis) {
    s=hjc::initialize(p); double e[6]={}; e[axis]=1.e-6;
    r=hjc::update(p,e,1.e-6,s);
    close(s.stress[axis],p.g*e[axis],1.e-6,"engineering shear convention");
    close(r.equivalent_stress,std::sqrt(3.)*p.g*e[axis],1.e-6,"von Mises shear invariant");
    e[axis]*=-1; hjc::update(p,e,1.e-6,s);
    close(s.stress[axis],0.,1.e-6,"elastic shear unloading");
  }
  // At zero pressure EFMIN bounds the denominator; it must not suppress damage.
  s=hjc::initialize(p); const double shear[6]={0,0,0,.01,0,0};
  r=hjc::update(p,shear,1.e-4,s);
  check(s.plastic_strain>0,"plastic return inactive");
  close(s.damage,s.plastic_strain/p.efmin,1.e-12,"low-pressure damage floor");
  p.d1=4.; s=hjc::initialize(p);
  const double initial=s.damage;
  close(initial,0.,0.,"initially undamaged material");
  for(int i=0;i<120;++i) {
    const double oldD=s.damage,oldEp=s.plastic_strain;
    r=hjc::update(p,shear,1.e-4,s);
    check(s.damage>=oldD && s.damage<=1 && s.plastic_strain>=oldEp,"irreversible history");
    check(std::isfinite(r.pressure) && r.acoustic_modulus>0,"finite result");
  }
  close(s.damage,1.,1.e-12,"damage saturation");
  const double compression[6]={-.01,-.01,-.01,.01,0,0};
  r=hjc::update(p,compression,1.e-4,s);
  check(r.pressure>0 && r.equivalent_stress>0,"damaged concrete retains confined strength");
  const double oldMup=s.plastic_volume;
  const double release[6]={.001,.001,.001,0,0,0};
  hjc::update(p,release,1.e-4,s);
  close(s.plastic_volume,oldMup,1.e-14,"plastic volume healed during unloading");
  // A known polynomial makes locking exactly mu=.12 (eta=.1, UL=1/55).
  p=parameters();p.pc=10.;p.muc=.01;p.pl=100.;p.mul=1./55.;
  p.k1=900.;p.k2=500.;p.k3=5000.;
  check(hjc::prepare(p),"continuous EOS parameters rejected");
  close(p.locking_compression,.12,1.e-14,"locking root");
  s=hjc::initialize(p);
  double e[6]={-std::log1p(.065)/3,-std::log1p(.065)/3,-std::log1p(.065)/3,0,0,0};
  r=hjc::update(p,e,1.e-4,s);
  close(r.pressure,55.,1.e-10,"crushing line midpoint");
  close(s.plastic_volume,1./110.,1.e-14,"crushing plastic volume midpoint");
  const double unload=(std::log1p(.065)-std::log1p(.04))/3;
  for(int i=0;i<3;++i) e[i]=unload;
  r=hjc::update(p,e,1.e-4,s);
  close(r.pressure,55.*(.04-1./110.)/(.065-1./110.),1.e-10,"unloading through plastic offset");
  s=hjc::initialize(p);
  for(int i=0;i<3;++i) e[i]=-std::log1p(.12)/3;
  r=hjc::update(p,e,1.e-4,s);
  close(r.pressure,100.,1.e-10,"continuous locking pressure");
  // Tensile cutoff has zero shear strength and still counts plastic flow.
  p=parameters();s=hjc::initialize(p);
  for(int i=0;i<3;++i) e[i]=-std::log1p(-2*p.tension/(p.pc/p.muc))/3;
  e[3]=.01;
  r=hjc::update(p,e,1.e-4,s);
  close(r.equivalent_stress,0.,1.e-7,"zero strength at tensile cutoff");
  check(s.plastic_strain>0 && s.damage>0,"tensile cutoff froze plastic flow");
  close(r.pressure,-p.tension*(1-s.damage),1.e-7,"updated tensile cutoff");
  // Very low rates cannot produce a negative radial scaling factor.
  p=parameters();p.c=1.; s=hjc::initialize(p);
  r=hjc::update(p,shear,1.e8,s);
  check(r.equivalent_stress>=0 && std::isfinite(r.equivalent_stress),"negative rate strength");
  close(r.equivalent_stress,p.a*p.fc,1.e-7,"sub-reference rate plateau");
}
}
#endif
