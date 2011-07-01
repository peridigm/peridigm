/*
 * calculators.h
 *
 *  Created on: Jun 20, 2011
 *      Author: jamitch
 */

#ifndef CALCULATORS_H_
#define CALCULATORS_H_

#include "utilities/Vector.h"
#include "mesh_input/quick_grid/QuickGrid.h"
#include "bond_volume_calculator.h"



namespace VOLUME_FRACTION {

using MATERIAL_EVALUATION::Bond_Volume_Calculator;
using UTILITIES::Vector3D;
using QUICKGRID::SpecRing2D;
using QUICKGRID::Spec1D;

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
	  axis((double[]){0,0,1.0})
	{
		/*
		 * This estimates the diagonal of a cell -- its not exact but probably good enough
		 */
		double ro=spec.getr0();
		diagonal=sqrt(DR*DR+(ro*D_THETA)*(ro*D_THETA)+DZ*DZ);
	}
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
};

} // namespace VOLUME_FRACTION

#endif /* CALCULATORS_H_ */
