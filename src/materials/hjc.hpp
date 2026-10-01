// SPDX-License-Identifier: BSD-3-Clause
// Holmquist-Johnson-Cook concrete; integration conventions are in examples/hjc.
#ifndef PERIDIGM_HJC_HPP
#define PERIDIGM_HJC_HPP

#include <cmath>
#ifdef __CUDACC__
#define HJC_HD __host__ __device__
#else
#define HJC_HD
#endif

namespace hjc {

struct Parameters {
  double g, a, b, c, n, fc, tension, rate0, efmin, sfmax;
  double pc, muc, pl, mul, d1, d2, k1, k2, k3;
  double locking_compression, crushing_slope;
};

// Stress ordering: xx, yy, zz, xy, yz, zx. All stresses are Cauchy stresses.
// Compression is positive for pressure and plastic_volume, negative for strain.
struct State {
  double stress[6];
  double volume_strain, plastic_volume, max_compression, density_measure;
  double damage, plastic_strain;
};

struct Result {
  // Consumed by the material adapter and material-point tests in other files.
  // cppcheck-suppress unusedStructMember
  double pressure, equivalent_stress, fracture_strain, acoustic_modulus;
};

inline bool valid(const Parameters& p) {
  const double values[] = {p.g,p.a,p.b,p.c,p.n,p.fc,p.tension,p.rate0,
    p.efmin,p.sfmax,p.pc,p.muc,p.pl,p.mul,p.d1,p.d2,p.k1,p.k2,p.k3};
  for(unsigned i=0;i<sizeof(values)/sizeof(values[0]);++i)
    if(!std::isfinite(values[i])) return false;
  return p.g>0 && p.a>=0 && p.b>=0 && p.c>=0 && p.n>0 && p.fc>0 &&
    p.tension>0 && p.rate0>0 && p.efmin>0 && p.sfmax>0 &&
    p.pc>0 && p.muc>0 && p.pl>p.pc && p.mul>0 && p.d1>0 && p.d2>=0 &&
    p.k1>0 && p.k2>=0 && p.k3>=0;
}

// Match the crushing line to the dense EOS at PL. UL is permanent compaction.
// Derived constants are prepared once, outside the material-point loop.
inline bool prepare(Parameters& p) {
  if(!valid(p)) return false;
  double lo=0.0,hi=p.pl/p.k1;
  for(int i=0;i<64;++i) {
    const double eta=0.5*(lo+hi);
    if(eta*(p.k1+eta*(p.k2+eta*p.k3))<p.pl) lo=eta;
    else hi=eta;
  }
  p.locking_compression=p.mul+(1.0+p.mul)*0.5*(lo+hi);
  if(p.locking_compression<=p.muc) return false;
  p.crushing_slope=(p.pl-p.pc)/(p.locking_compression-p.muc);
  return std::isfinite(p.crushing_slope) && std::isfinite(p.pc/p.muc) &&
    p.crushing_slope<fmin(p.pc/p.muc,p.pl/(p.locking_compression-p.mul));
}

HJC_HD inline State initialize(const Parameters&) {
  State s={};
  return s;
}

// deps contains engineering shear increments. The caller supplies an objective
// frame and manages initialization, rotations, timestep acceptance and erosion.
HJC_HD inline Result update(const Parameters& p, const double deps[6],
                            double dt, State& s) {
  const double third=(deps[0]+deps[1]+deps[2])/3.0;
  double dev[6];
  for(int i=0;i<3;++i) dev[i]=deps[i]-third;
  for(int i=3;i<6;++i) dev[i]=0.5*deps[i];
  const double norm=dev[0]*dev[0]+dev[1]*dev[1]+dev[2]*dev[2]+
    2.0*(dev[3]*dev[3]+dev[4]*dev[4]+dev[5]*dev[5]);
  const double rate=dt>0 ? sqrt(2.0/3.0)*sqrt(norm)/dt : 0.0;
  const double oldp=-(s.stress[0]+s.stress[1]+s.stress[2])/3.0;
  for(int i=0;i<3;++i) s.stress[i]+=oldp+2.0*p.g*dev[i];
  for(int i=3;i<6;++i) s.stress[i]+=p.g*deps[i];
  const double trial=sqrt(1.5*(s.stress[0]*s.stress[0]+
    s.stress[1]*s.stress[1]+s.stress[2]*s.stress[2])+
    3.0*(s.stress[3]*s.stress[3]+s.stress[4]*s.stress[4]+
    s.stress[5]*s.stress[5]));

  s.volume_strain+=deps[0]+deps[1]+deps[2];
  const double mu=expm1(-s.volume_strain);
  const double k0=p.pc/p.muc;
  double pressure;
  s.max_compression=fmax(s.max_compression,mu);
  s.density_measure=fmin(1.0,fmax(0.0,
    (s.max_compression-p.muc)/(p.locking_compression-p.muc)));
  const double volume=p.mul*s.density_measure;
  const double dmup=fmax(0.0,volume-s.plastic_volume);
  s.plastic_volume=volume;
  if(s.max_compression>=p.locking_compression) {
    const double m=(mu-p.mul)/(1.0+p.mul);
    pressure=p.k1*m;
    if(m>0.0) pressure+=m*m*(p.k2+p.k3*m);
  } else if(s.max_compression>p.muc) {
    const double peakPressure=p.pc+p.crushing_slope*(s.max_compression-p.muc);
    pressure=peakPressure/(s.max_compression-volume)*(mu-volume);
  } else {
    pressure=k0*mu;
  }
  pressure=fmax(pressure,-p.tension*(1.0-s.damage));

  double strength=pressure<0.0 ?
    fmax(0.0,p.a*(1.0-s.damage+pressure/p.tension)) :
    p.a*(1.0-s.damage)+p.b*pow(pressure/p.fc,p.n);
  strength=fmax(0.0,fmin(p.sfmax,
    strength*(1.0+p.c*log(fmax(1.0,rate/p.rate0)))));
  const double normalized_trial=trial/p.fc;
  const double scale=normalized_trial<strength ? 1.0 :
    strength/fmax(normalized_trial,1.e-12);
  for(int i=0;i<6;++i) s.stress[i]*=scale;
  const double dep=(1.0-scale)*trial/(3.0*p.g);
  s.plastic_strain+=dep;
  const double ef=fmax(p.efmin,p.d1*pow(fmax((pressure+p.tension)/p.fc,0.0),p.d2));
  s.damage=fmin(1.0,s.damage+(dep+dmup)/ef);
  pressure=fmax(pressure,-p.tension*(1.0-s.damage));
  for(int i=0;i<3;++i) s.stress[i]-=pressure;
  // Bound wave stiffness using the largest compression reached so far.
  const double eta=fmax(0.0,(s.max_compression-p.mul)/(1.0+p.mul));
  const double dense=(p.k1+2.0*p.k2*eta+3.0*p.k3*eta*eta)/(1.0+p.mul);
  Result r={pressure,trial*scale,ef,
    (1.0+fmax(0.0,s.max_compression))*fmax(k0,dense)+4.0*p.g/3.0};
  return r;
}

} // namespace hjc
#undef HJC_HD
#endif
