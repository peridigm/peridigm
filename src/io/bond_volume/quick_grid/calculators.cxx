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

#include "calculators.h"
#include "utilities/Array.h"

#include <iostream>
#include <stdexcept>
#include <string>

#include <Teuchos_ParameterList.hpp>
#ifdef USE_YAML
  #include <Teuchos_YamlParameterListCoreHelpers.hpp>
#endif

namespace BOND_VOLUME {

namespace QUICKGRID {

using UTILITIES::Array;

  std::shared_ptr<Bond_Volume_Calculator> get_Bond_Volume_Calculator(const std::string& yaml_file_name) {

   Teuchos::ParameterList params;
   Teuchos::Ptr<Teuchos::ParameterList> params_ptr(&params);
   if (!::QUICKGRID::file_exists(yaml_file_name)) {
     std::string msg = "\n**** Error, failed to open file " + yaml_file_name;
     msg += " (unit tests may need to be executed from within the directory where they reside).";
     std::cout << "\n" << msg << std::endl;;
     throw std::runtime_error(msg);
   }
#ifdef USE_YAML
   Teuchos::updateParametersFromYamlFile(yaml_file_name, params_ptr);
#else
   std::string msg = "\n\n**** Error, get_Bond_Volume_Calculator() requires YAML.";
   std::cout << msg << std::endl;;
   throw std::runtime_error(msg);
#endif

   std::shared_ptr<Bond_Volume_Calculator> c;

   Teuchos::ParameterList disc_params = params.sublist("Discretization");
   double horizon = disc_params.get<double>("Horizon");
   std::string type = disc_params.get<std::string>("Type");

   if(type == "QuickGrid.TensorProduct3DMeshGenerator"){
    Teuchos::ParameterList mesh_params = params.sublist("QuickGrid").sublist("TensorProduct3DMeshGenerator");
   	double xStart = mesh_params.get<double>("X Origin");
   	double yStart = mesh_params.get<double>("Y Origin");
   	double zStart = mesh_params.get<double>("Z Origin");
   	double xLength = mesh_params.get<double>("X Length");
   	double yLength = mesh_params.get<double>("Y Length");
   	double zLength = mesh_params.get<double>("Z Length");
   	const int nx = mesh_params.get<int>("Number Points X");
   	const int ny = mesh_params.get<int>("Number Points Y");
   	const int nz = mesh_params.get<int>("Number Points Z");
		const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
		const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
		const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
		c = std::shared_ptr<VolumeFractionCalculator>(new VolumeFractionCalculator(xSpec,ySpec,zSpec,horizon));
	} else {
		std::string s;
		s = "Error-->BOND_VOLUME::QUICKGRID::get_Bond_Volume_Calculator()\n";
		s += "\tOnly Reader for Discretization.Type==QuickGrid.TensorProduct3DMeshGenerator is implemented\n";
		s += "\tCome back soon for the other type(s):)\n";
		throw std::runtime_error(s);
	}

	return c;

}


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
		int numNeigh = *neighPtr; neighPtr++;

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
	int numNeigh = *neighPtr; neighPtr++;
	for(int n=0;n<numNeigh;n++,neighPtr++,bond_volumes++){
		int localId = *neighPtr;
		const double *Q = &xOverlap[3*localId];
		*bond_volumes = c->operator ()(P,Q);
	}




}


} // namespace QUICKGRID

} // namespace VOLUME_FRACTION
