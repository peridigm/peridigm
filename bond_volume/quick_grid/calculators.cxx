/*
 * calculators.cxx
 *
 *  Created on: Jun 20, 2011
 *      Author: jamitch
 */
#include "calculators.h"
#include "utilities/Array.h"
#include <stdexcept>



namespace BOND_VOLUME {

namespace QUICKGRID {

using UTILITIES::Array;

double RingVolumeFractionCalculator::cellVolume(const double* q) const {
	Vector3D Q(*q,*(q+1),*(q+2));
	Array<double> r, theta, z;
	double dr, d_theta, dz;
	{
		// Calculate centroid of cell
		/*
		 * rC -- radial coordinate of cell centroid
		 * thetaC -- theta coordinate of cell centroid
		 * zC -- z coordinate of cell centroid
		 */
		UTILITIES::Minus minus;
		UTILITIES::Dot dot;
		// qc=q-c
		Vector3D qc = minus(Q,c);
		// height of point q relative to c along axis of cylinder
		double h = dot(qc,axis);
		// vector of height along axis
		Vector3D qZ(axis); qZ*=h;
		// projection of qc onto r-theta plane
		Vector3D qR = minus(qc,qZ);
		double rC = qR.norm();
		double x = qR[0];
		double y = qR[1];
		double zC = c[2]+h;
		// angular centroid of point c
		double thetaC = atan2(y, x);
		//		double dv = rC * dr * d_theta * dz;
		Spec1D rSpec(nR,rC-DR/2,DR);
		Spec1D thetaSpec(nTheta,thetaC-D_THETA/2,D_THETA);
		Spec1D zSpec(nZ,zC-DZ/2,DZ);
		r = getDiscretization(rSpec);
		theta = getDiscretization(thetaSpec);
		z = getDiscretization(zSpec);
		dr = rSpec.getCellSize();
		d_theta = thetaSpec.getCellSize();
		dz = zSpec.getCellSize();

	}
	/*
	 * check each point 'q' of discretized cell
	 */
	double volume=0;
	for(size_t i=0;i<nR;i++){
		double dV = r[i] * dr * d_theta * dz;
		for(size_t j=0;j<nTheta;j++){
			for(size_t k=0;k<nZ;k++)
				volume +=dV;
		}

	}
	return volume;
}


double RingVolumeFractionCalculator::operator() (const double* pCenter, const double* qNeigh) const {

	const Vector3D P(*pCenter,*(pCenter+1),*(pCenter+2));
	const Vector3D Q(*qNeigh,*(qNeigh+1),*(qNeigh+2));
	Array<double> r, theta, z;
	double dr, d_theta, dz;
	{
		// Calculate centroid of cell
		/*
		 * rC -- radial coordinate of cell centroid
		 * thetaC -- theta coordinate of cell centroid
		 * zC -- z coordinate of cell centroid
		 */
		UTILITIES::Minus minus;
		UTILITIES::Dot dot;
		// qc=q-c
		Vector3D qc = minus(Q,c);
		// height of point q relative to c along axis of cylinder
		double h = dot(qc,axis);
		// vector of height along axis
		Vector3D qZ(axis); qZ*=h;
		// projection of qc onto r-theta plane
		Vector3D qR = minus(qc,qZ);
		double rC = qR.norm();
		double x = qR[0];
		double y = qR[1];
		double zC = c[2]+h;
		// angular centroid of point c
		double thetaC = atan2(y, x);
		//		double dv = rC * dr * d_theta * dz;
		Spec1D rSpec(nR,rC-DR/2,DR);
		Spec1D thetaSpec(nTheta,thetaC-D_THETA/2,D_THETA);
		Spec1D zSpec(nZ,zC-DZ/2,DZ);
		r = getDiscretization(rSpec);
		theta = getDiscretization(thetaSpec);
		z = getDiscretization(zSpec);
		dr = rSpec.getCellSize();
		d_theta = thetaSpec.getCellSize();
		dz = zSpec.getCellSize();

	}
	/*
	 * check each point 'q' of discretized cell
	 */
	double volume=0;
	Vector3D q;
	for(size_t i=0;i<nR;i++){
		double dV = r[i] * dr * d_theta * dz;
		for(size_t j=0;j<nTheta;j++){
			double x = r[i]*cos(theta[j]);
			double y = r[i]*sin(theta[j]);
			for(size_t k=0;k<nZ;k++){
				q[0]=x;
				q[1]=y;
				q[2]=z[k];
				volume += (comparator(P,q) ? dV : 0.0);
			}
		}

	}
	return volume;
}

double VolumeFractionCalculator::operator() (const double* pCenter, const double* qNeigh) const {

	const Vector3D P(*pCenter,*(pCenter+1),*(pCenter+2));
	const Vector3D Q(*qNeigh,*(qNeigh+1),*(qNeigh+2));
	Spec1D xSpec(nX,Q[0]-DX/2,DX);
	Spec1D ySpec(nY,Q[1]-DY/2,DY);
	Spec1D zSpec(nZ,Q[2]-DZ/2,DZ);
	double dx(xSpec.getCellSize()), dy(ySpec.getCellSize()), dz(zSpec.getCellSize());
	double dV=dx*dy*dz;
	Array<double> xArray, yArray, zArray;
	xArray = getDiscretization(xSpec);
	yArray = getDiscretization(ySpec);
	zArray = getDiscretization(zSpec);

	/*
	 * check each point 'q' of discretized cell
	 */
	double volume=0;
	Vector3D q;
	for(size_t i=0;i<nX;i++){
		double x = xArray[i];
		for(size_t j=0;j<nY;j++){
			double y = yArray[j];
			for(size_t k=0;k<nZ;k++){
				q[0]=x;
				q[1]=y;
				q[2]=zArray[k];
				volume += (comparator(P,q) ? dV : 0.0);
			}
		}

	}
	return volume;
}


void compute_bond_volume
(
		size_t num_owned_points,
		const int*  localNeighborList,
		const double* xOverlap,
		double *bond_volumes,
		const Bond_Volume_Calculator *c
) {
	const int *neighPtr = localNeighborList;
	const double *xOwned = xOverlap;
	for(size_t p=0;p<num_owned_points;p++, xOwned +=3){
		size_t numNeigh = *neighPtr; neighPtr++;

		const double *P = xOwned;

		/*
		 * Loop over neighborhood of point P and compute
		 * bond volume associated with each point Q
		 */
		for(int n=0;n<numNeigh;n++,neighPtr++,bond_volumes++){
			int localId = *neighPtr;
			const double *Q = &xOverlap[3*localId];
			*bond_volumes = c->operator ()(P,Q);
		}
	}

}

void compute_bond_volume
(
		const double *X,
		const int*  localNeighborList,
		const double* xOverlap,
		double *bond_volumes,
		const Bond_Volume_Calculator *c
) {

	const int *neighPtr = localNeighborList;
	const double *P = X;
	size_t numNeigh = *neighPtr; neighPtr++;
	for(int n=0;n<numNeigh;n++,neighPtr++,bond_volumes++){
		int localId = *neighPtr;
		const double *Q = &xOverlap[3*localId];
		*bond_volumes = c->operator ()(P,Q);
	}




}


} // namespace QUICKGRID

} // namespace VOLUME_FRACTION
