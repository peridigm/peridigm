//! \file elastic_kokkos.h

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
#ifndef ELASTIC_KOKKOS_H
#define ELASTIC_KOKKOS_H


// TODO: (DZT) Move most of this to a typedefs file

#ifdef _OPENMP
  #include <Kokkos_OpenMP.hpp>
#else
  #include <Kokkos_Threads.hpp>
#endif

#ifdef KOKKOS_HAVE_CUDA
  #include <Kokkos_Cuda.hpp>
  #include <cuda.h>
  #include <cuda_runtime.h>
  typedef Kokkos::Cuda device_type;
#else
  #ifdef _OPENMP
    typedef Kokkos::OpenMP device_type;
  #else
    typedef Kokkos::Threads device_type;
  #endif

#endif


#include <Kokkos_View.hpp>
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Timer.hpp>

/* Define types used throughout the code */

//Position arrays
typedef Kokkos::View<double*[3], Kokkos::LayoutRight, device_type>                                   t_x_array ;
typedef t_x_array::HostMirror                                                                        t_x_array_host ;
typedef Kokkos::View<const double*[3], Kokkos::LayoutRight, device_type>                             t_x_array_const ;
typedef Kokkos::View<const double*[3], Kokkos::LayoutRight, device_type, Kokkos::MemoryRandomAccess >  t_x_array_randomread ;

//Force array
//typedef Kokkos::View<double*[3],  device_type>                                                       t_f_array ;

// TODO REMove unneccesary typdefs


//Neighborlist
typedef Kokkos::View<int**, device_type >                                                            t_neighbors ;
typedef Kokkos::View<const int**, device_type >                                                      t_neighbors_const ;
typedef Kokkos::View<int*, device_type, Kokkos::MemoryUnmanaged >                                    t_neighbors_sub ;
typedef Kokkos::View<const int*, device_type, Kokkos::MemoryUnmanaged >                              t_neighbors_const_sub ;

//1d int array
typedef Kokkos::View<int*, device_type >                                                             t_int_1d ;
typedef t_int_1d::HostMirror                                                                         t_int_1d_host ;
typedef Kokkos::View<const int*, device_type >                                                       t_int_1d_const ;
typedef Kokkos::View<int*, device_type , Kokkos::MemoryUnmanaged>                                    t_int_1d_um ;
typedef Kokkos::View<const int* , device_type , Kokkos::MemoryUnmanaged>                             t_int_1d_const_um ;

//2d int array
typedef Kokkos::View<int**, Kokkos::LayoutRight, device_type >                                       t_int_2d ;
typedef t_int_2d::HostMirror                                                                         t_int_2d_host ;

//Scalar ints
typedef Kokkos::View<int[1], Kokkos::LayoutLeft, device_type>                                        t_int_scalar ;
typedef t_int_scalar::HostMirror                                                                     t_int_scalar_host ;

typedef Kokkos::View<double*, device_type>                                                           t_double_1d ;
typedef t_double_1d::HostMirror                                                                      t_double_1d_host ;
typedef Kokkos::View<const double*, device_type, Kokkos::MemoryRandomAccess >                        t_double_1d_randomread ;


namespace MATERIAL_EVALUATION {

struct System {
  int nlocal;
  int noverlap;
  //double horizon;
  double thermal_coef;
  double K;
  double MU;

  // TODO change t_double to be generic
  t_double_1d d_x;
  t_double_1d d_y;
  t_double_1d d_theta;
  t_double_1d d_delta_t;
  t_double_1d d_vol;
  t_double_1d d_m;
  t_double_1d_host h_x;
  t_double_1d_host h_y;
  t_double_1d_host h_theta;
  t_double_1d_host h_delta_t;
  t_double_1d_host h_vol;
  t_double_1d_host h_m;

  t_double_1d d_f;
  t_double_1d_host h_f;

  t_neighbors neighbors;
  t_neighbors bond_damage;
  t_int_1d numneigh;
};



int create_system(System &system, int nx, int ny, int nz, double rho);
void neigh_setup(System &system);
void neigh_build(System &system);
void force(System &system);


template<typename ScalarT>
void computeInternalForceLinearElasticKokkos
(
    const double* xOverlapPtr,
    const ScalarT* yOverlapPtr,
    const double* mOwned,
    const double* volumeOverlapPtr,
    const ScalarT* dilatationOwned,
    const double* bondDamage,
    const double* dsfOwned,
    ScalarT* fInternalOverlapPtr,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double SHEAR_MODULUS,
        double horizon,
        double thermalExpansionCoefficient = 0,
        const double* deltaTemperature = 0

);


} // MATERIAL_EVALUATION



#endif // ELASTIC_KOKKOS_H
