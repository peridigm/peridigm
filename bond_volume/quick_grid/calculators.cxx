/*
 * calculators.cxx
 *
 *  Created on: Jun 20, 2011
 *      Author: jamitch
 */
#include "calculators.h"
#include "utilities/Array.h"
#include <stdexcept>



namespace VOLUME_FRACTION {

using UTILITIES::Array;

double RingVolumeFractionCalculator::cellVolume(const double* q) const {
	Vector3D Q;
        Q[0] = q[0]; Q[1] = q[1]; Q[2] = q[2];
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

	Vector3D P,Q;
        P[0] = pCenter[0]; P[1] = pCenter[1]; P[2] = pCenter[2];
        Q[0] = qNeigh[0]; Q[1] = qNeigh[1]; Q[2] = qNeigh[2];
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

	//const Vector3D P((double[3]){*pCenter,*(pCenter+1),*(pCenter+2)});
	//const Vector3D Q((double[3]){*qNeigh,*(qNeigh+1),*(qNeigh+2)});
	Vector3D P,Q;
        P[0] = pCenter[0]; P[1] = pCenter[1]; P[2] = pCenter[2];
        Q[0] = qNeigh[0]; Q[1] = qNeigh[1]; Q[2] = qNeigh[2];
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

} // namespace VOLUME_FRACTION
