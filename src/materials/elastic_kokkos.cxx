//! \file elastic_kokkos.cxx

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
#include "elastic_kokkos.h"
#include <cmath>
#include <Sacado.hpp>
#include "material_utilities.h"

/* Define values which set the max number of registers used for the Force Kernel
 * Its 32 * 2048 / (KOKKOS_CUDA_MAX_THREADS * KOKKOS_CUDA_MIN_BLOCKS)
 * Have to be set before including Kokkos header files.
 */

#define KOKKOS_CUDA_MAX_THREADS 512
#define KOKKOS_CUDA_MIN_BLOCKS 3

#include <Kokkos_Atomic.hpp>
#include <cstdio>
#include <cstdlib>

#define SMALL 1.0e-6
#define FACTOR 0.999
/* initialize atoms on fcc lattice in parallel fashion */

#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a<b?a:b)

namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void computeInternalForceLinearElasticKokkos
(
    const double* xOverlap,
    const ScalarT* yOverlap,
    const double* mOwned,
    const double* volumeOverlap,
    const ScalarT* dilatationOwned,
    const double* bondDamage,
    const double* dsfOwned,
    ScalarT* fInternalOverlap,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double SHEAR_MODULUS,
        double horizon,
        double thermalExpansionCoefficient,
        const double* deltaTemperature
)
{
   /* Thread numbers for Host */

   // TODO bring this in from input/choose an appropriate thread count
   int num_threads = 1;
   int teams = 1;
   int device = 0; // Default device for GPU runs

 #ifdef KOKKOS_HAVE_CUDA
   Kokkos::Cuda::host_mirror_device_type::initialize(teams*num_threads);
   Kokkos::Cuda::SelectDevice select_device(device);
   Kokkos::Cuda::initialize(select_device);
 #else
   #ifdef _OPENMP
   Kokkos::OpenMP::initialize(teams*num_threads);
   #else
   Kokkos::Threads::initialize(teams*num_threads);
   #endif
 #endif

   const int dim = 3;

   System system;

   system.nlocal       = numOwnedPoints;
   system.thermal_coef = thermalExpansionCoefficient;
   system.numneigh     = t_int_1d("numneigh",system.nlocal);
   system.K            = BULK_MODULUS;
   system.MU           = SHEAR_MODULUS;
   system.d_theta      = t_double_1d("theta",system.nlocal);
   system.d_delta_t    = t_double_1d("delta_t",system.nlocal);
   system.d_m          = t_double_1d("m",system.nlocal);
   system.h_theta      = Kokkos::create_mirror_view(system.d_theta);
   system.h_delta_t    = Kokkos::create_mirror_view(system.d_delta_t);
   system.h_m          = Kokkos::create_mirror_view(system.d_m);

   // copy the data from the overlap vecs to the kokkos views
   // TODO (DZT) api to the epetra vectors with the view to skip this step (check with kokkos guys in this)
   const int *neighPtrHost  = localNeighborList;

   int maxNumNeigh = 0;
   for(int p=0;p<numOwnedPoints;p++){
     int numNeigh = *neighPtrHost; neighPtrHost++;
     system.numneigh[p] = numNeigh;
     if(numNeigh>maxNumNeigh)
       maxNumNeigh = numNeigh;
     neighPtrHost+= numNeigh;
   }

   const double *mHost      = mOwned;
   const ScalarT *thetaHost = dilatationOwned;
   const double *deltaTHost = deltaTemperature;

   system.neighbors = t_neighbors("neighbors",system.nlocal,maxNumNeigh);
   system.bond_damage = t_neighbors("bond_damage",system.nlocal,maxNumNeigh);
   int maxNeighIndex = 0;
   const int *neighPtrList = localNeighborList;
   const double *damagePtr = bondDamage;
   for(int p=0;p<numOwnedPoints;p++,thetaHost++,deltaTHost++,mHost++){
     system.h_m[p] = *mHost;
     system.h_theta[p] = *thetaHost;
     system.h_delta_t[p] = (deltaTemperature) ? *deltaTHost : 0.0;

     int numNeigh = *neighPtrList; neighPtrList++;
     for(int n=0;n<numNeigh;n++,neighPtrList++,damagePtr++){
       system.neighbors(p,n) = *neighPtrList;
       system.bond_damage(p,n) = *damagePtr;
       if(*neighPtrList > maxNeighIndex)
         maxNeighIndex = *neighPtrList;
     }
   }
   system.noverlap     = maxNeighIndex + 1;
   system.d_x          = t_double_1d("X",system.noverlap*dim);
   system.d_y          = t_double_1d("Y",system.noverlap*dim);
   system.d_vol        = t_double_1d("vol",system.noverlap);
   system.d_f          = t_double_1d("F",system.noverlap*dim,0.0);
   system.h_x          = Kokkos::create_mirror_view(system.d_x);
   system.h_y          = Kokkos::create_mirror_view(system.d_y);
   system.h_vol        = Kokkos::create_mirror_view(system.d_vol);
   system.h_f          = Kokkos::create_mirror_view(system.d_f);

   // copy the data from the overlap vecs to the kokkos views
   // TODO (DZT) api to the epetra vectors with the view to skip this step (check with kokkos guys on this)
   const double *xHost      = xOverlap;
   const double *yHost      = yOverlap;
   const int *neighPtrHost2 = localNeighborList;
   const double *vHost      = volumeOverlap;

   for(int p=0;p<system.noverlap;p++,xHost +=dim,yHost +=dim, vHost++){
     const double *X = xHost;
     const double *Y = yHost;
     system.h_x[p*dim+0] = X[0];
     system.h_x[p*dim+1] = X[1];
     system.h_x[p*dim+2] = X[2];
     system.h_y[p*dim+0] = Y[0];
     system.h_y[p*dim+1] = Y[1];
     system.h_y[p*dim+2] = Y[2];
     system.h_vol[p] = *vHost;
   }

   Kokkos::deep_copy(system.d_x,system.h_x); // TODO (DZT) Don't need to deep copy this each time
   Kokkos::deep_copy(system.d_y,system.h_y);
   Kokkos::deep_copy(system.d_vol,system.h_vol); // TODO (DZT) Don't need to deep copy this each time
   Kokkos::deep_copy(system.d_theta,system.h_theta);
   Kokkos::deep_copy(system.d_delta_t,system.h_delta_t);
   Kokkos::deep_copy(system.d_m,system.h_m); // TODO (DZT) Don't need to deep copy this each time

   force(system);

   device_type::finalize();

  // copy the force back over:
  for(int p=0;p<system.noverlap;p++){
    for(int i=0;i<dim;++i)
      fInternalOverlap[p*dim+i] = system.h_f[p*dim+i];
  }
}
/** Explicit template instantiation for double. */
template void computeInternalForceLinearElasticKokkos<double>
(
    const double* xOverlap,
    const double* yOverlap,
    const double* mOwned,
    const double* volumeOverlap,
    const double* dilatationOwned,
    const double* bondDamage,
    const double* dsfOwned,
    double* fInternalOverlap,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double SHEAR_MODULUS,
        double horizon,
        double thermalExpansionCoefficient,
        const double* deltaTemperature
 );


struct ForceFunctor {

  typedef t_x_array::device_type device_type; //Device Type for running the kernel

  t_double_1d_randomread x;      //positions
  t_double_1d_randomread y;      //current positions
  t_double_1d_randomread vol;    //cell volume
  t_double_1d_randomread m;      //weighted volume
  t_double_1d_randomread theta;  //dilatation
  t_double_1d_randomread delta_t;//temperature change
  t_int_1d_const numneigh;       //number of neighbors per cell
  t_neighbors_const neighbors;   //neighborlist
  t_neighbors_const bond_damage; //damage
  System sys;
  double K;
  double MU;
  double thermal_coef;
  //double horizon;

  ForceFunctor(System & s) {
    sys = s;
    // TODO (DZT) convert all these to point to system
    x = s.d_x;
    y = s.d_y;
    theta = s.d_theta;
    delta_t = s.d_delta_t;
    vol = s.d_vol;
    m = s.d_m;
    numneigh = s.numneigh;
    neighbors = s.neighbors;
    bond_damage = s.bond_damage;
    K = s.K;
    MU = s.MU;
    thermal_coef = s.thermal_coef;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i) const {
    force(i);
  }

  KOKKOS_INLINE_FUNCTION
  void force(const int &i) const
  {
    const int numneighs = numneigh[i];
    const double Xtmp = x[i*3+0];
    const double Ytmp = x[i*3+1];
    const double Ztmp = x[i*3+2];
    const double xtmp = y[i*3+0];
    const double ytmp = y[i*3+1];
    const double ztmp = y[i*3+2];
    const double cellVolumeSelf = vol[i];
    const double deltaTSelf = delta_t[i];
    const double thetaSelf = theta[i];
    const double MSelf = m[i];
    const double alphaSelf = 15.0*MU/MSelf;

    double cellVolumeNeigh, delX, delY, delZ, delx, dely, delz, zeta;
    double omega, dY, t, fx, fy, fz, e, c1;

    for(int k=0; k < numneighs; k++) {
      const int j = neighbors(i, k);
      cellVolumeNeigh = vol[j];

      delX = x[j*3+0] - Xtmp;
      delY = x[j*3+1] - Ytmp;
      delZ = x[j*3+2] - Ztmp;
      delx = y[j*3+0] - xtmp;
      dely = y[j*3+1] - ytmp;
      delz = y[j*3+2] - ztmp;
      zeta = sqrt(delX*delX+delY*delY+delZ*delZ);
      dY = sqrt(delx*delx+dely*dely+delz*delz);
      e = dY - zeta - thermal_coef*deltaTSelf*zeta;
      omega = 1.0; // TODO add influence function
      //omega = scalarInfluenceFunction(zeta,horizon);
      c1 = omega*thetaSelf*(3.0*K/MSelf-alphaSelf/3.0);
      t = (1.0-bond_damage(i,k))*(c1 * zeta + (1.0-bond_damage(i,k)) * omega * alphaSelf * e);
      fx = t * delx / dY;
      fy = t * dely / dY;
      fz = t * delz / dY;

      sys.d_f[i*3+0] += fx * cellVolumeNeigh;
      sys.d_f[i*3+1] += fy * cellVolumeNeigh;
      sys.d_f[i*3+2] += fz * cellVolumeNeigh;
      sys.d_f[j*3+0] -= fx * cellVolumeSelf;
      sys.d_f[j*3+1] -= fy * cellVolumeSelf;
      sys.d_f[j*3+2] -= fz * cellVolumeSelf;
    }
  }
};

/* Calling function */

void force(System &s){

  ForceFunctor f(s);
  Kokkos::parallel_for(s.nlocal,f);
  device_type::fence();
}

} // MATERIAL_EVALUATION
