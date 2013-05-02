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

#ifndef CALCULATORS_H_
#define CALCULATORS_H_

#include "utilities/Vector3D.h"
#include "mesh_input/quick_grid/QuickGrid.h"
#include "bond_volume/bond_volume_calculator.h"



namespace BOND_VOLUME {

namespace QUICKGRID {

using MATERIAL_EVALUATION::Bond_Volume_Calculator;
using UTILITIES::Vector3D;
using ::QUICKGRID::SpecRing2D;
using ::QUICKGRID::Spec1D;


struct RingVolumeFractionCalculator : public Bond_Volume_Calculator  {

private:
	UTILITIES::InsideSphere comparator;
	const double DR, D_THETA, DZ;
	const size_t nR, nTheta, nZ;
	const Vector3D c, axis;
	double diagonal;

public:
	/*
	 * NOTE BIG ASSUMPTION that is EASY TO FIX
	 * CYLINDER AXIS is along z-coordinate
	 *
	 * May need to check assumptions about center in spec
	 *     if non-trivial origin is used
	 */
	RingVolumeFractionCalculator(const SpecRing2D& spec, const Spec1D& axisSpec, double horizon)
	: Bond_Volume_Calculator(),
	  comparator(horizon),
	  DR(spec.getRaySpec().getCellSize()),
	  D_THETA(spec.getRingSpec().getCellSize()),
	  DZ(axisSpec.getCellSize()),
	  nR(16),
	  nTheta(16),
	  nZ(16),
	  c(spec.getCenter()),
	  axis(0.0,0.0,1.0)
	{
		/*
		 * This estimates the diagonal of a cell -- its not exact but probably good enough
		 */
		double ro=spec.getr0();
		diagonal=sqrt(DR*DR+(ro*D_THETA)*(ro*D_THETA)+DZ*DZ);
	}
	virtual ~RingVolumeFractionCalculator(){}
	double get_cell_diagonal() const { return diagonal; }
	/**
	 * Q is neighbor of P
	 * This function computes volume contribution of Q to P neighborhood
	 */
	double operator() (const double* p, const double* q) const;
	/*
	 * Compute the volume of cell Q using QUICKGRID quadrature
	 */
	double cellVolume(const double* q) const;

	double get_horizon() const { return comparator.get_radius(); }
};


struct VolumeFractionCalculator : public Bond_Volume_Calculator  {

private:
	UTILITIES::InsideSphere comparator;
	const double DX, DY, DZ;
	const size_t nX, nY, nZ;
	double diagonal;

public:
	VolumeFractionCalculator(const Spec1D& xSpec, const Spec1D& ySpec, const Spec1D& zSpec, double horizon)
	: Bond_Volume_Calculator(),
	  comparator(horizon),
	  DX(xSpec.getCellSize()),
	  DY(ySpec.getCellSize()),
	  DZ(zSpec.getCellSize()),
	  nX(16),
	  nY(16),
	  nZ(16)
	{
		/*
		 * This estimates the diagonal of a cell
		 */
		diagonal=sqrt(DX*DX+DY*DY+DZ*DZ);
	}
	virtual ~VolumeFractionCalculator(){}
	double get_cell_diagonal() const { return diagonal; }
	/**
	 * Q is neighbor of P
	 * This function computes volume contribution of Q to P neighborhood
	 */
	double operator() (const double* p, const double* q) const;
	/*
	 * Compute the volume of cell Q using QUICKGRID quadrature
	 */
	double cellVolume(const double* q) const { return DX * DY * DZ; }

	double get_horizon() const { return comparator.get_radius(); }
};


/**
 * Call this function to compute bond volumes.  Those
 * are volumes relating to a center point P and its neighbors
 * points Q.  At the outer reaches of a neighborhood horizon,
 * cells Q will typically only contribute a fraction of their
 * cell volume.
 * Input:
 *  num_owned_points -- number of points owned by this processor
 *  localNeighborhoodList -- this is the usually neighborhood list
 *          that has the number of neighbors associated with each
 *          owned point included at the start of each list
 *  xOverlap -- mesh coordinates which includes overlap/ghosted points
 *  Bond_Volume_Calculator: this is the operator that evaluates each bond
 *          volume
 *  Output:
 *   bond_volumes -- this will be calculated but incoming pointer must
 *          be pre-allocated.  Data layout is identical to the neighborhood list
 *   			except that it does not include an extra value at the start (num neighbors)
 *
 */
void compute_bond_volume
(
		size_t num_owned_points,
		const int*  localNeighborList,
		const double* xOverlap,
		double *bond_volumes,
		const Bond_Volume_Calculator *c
);

/**
 * Call this function on a single point 'X'
 * NOTE: neighPtr to should point to 'numNeigh' for 'X'
 * and thus describe the neighborhood list as usual
 */
void compute_bond_volume
(
		const double *X,
		const int*  localNeighborList,
		const double* xOverlap,
		double *bond_volumes,
		const Bond_Volume_Calculator *c
);

/**
 * prototype
 */
 std::tr1::shared_ptr<Bond_Volume_Calculator> get_Bond_Volume_Calculator(const std::string& json_filename);

} // namespace QUICKGRID

} // namespace VOLUME_FRACTION

#endif /* CALCULATORS_H_ */
