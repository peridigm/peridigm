//! \file rkpmshapefunction.cxx

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
//#include <iostream>
#include <cmath>
#include <Sacado.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialQRDenseSolver.hpp>
#include "rkpmshapefunction.h"

namespace RKPM_EVALUATION {

void computeShapeFunction
(
		const double* xOverlap,
		const double* bondDamage,
		int myNumPoints,
		const int*  localNeighborList,
		double* OwnRKPMShapeFunctionValues,
		double* BondRKPMShapeFunctionValues,
        double RKPMSupport,
        int RKPMBasisOrder,
        double PHI(double)
)
{

	/*
	 * Compute RKPM Shape Function
	 */
    ///*
    double xi[3], xn[3];
    int numNeighbors, neighborIndex, neighborListIndex(0), rkpmFunctionIndex(0);
    for(int p=0 ; p<myNumPoints ; p++)
    {
      xi[0] = xOverlap[3*p];
      xi[1] = xOverlap[3*p+1];
      xi[2] = xOverlap[3*p+2];
      numNeighbors = localNeighborList[neighborListIndex++];


      //int RKDimension = RKPMBasisOrder+1; //1D
      //int RKDimension = (RKPMBasisOrder+1)*(RKPMBasisOrder+2)/2; //2D
      int RKDimension = (RKPMBasisOrder+1)*(RKPMBasisOrder+2)*(RKPMBasisOrder+3)/6; //3D

      Teuchos ::SerialDenseMatrix<int, double> M(RKDimension, RKDimension);
      Teuchos ::SerialDenseMatrix<int, double> b(RKDimension, 1);
      Teuchos ::SerialDenseMatrix<int, double> hRK(RKDimension, 1);

      double phi;
      double distanceKernel = 0.0;

      M.scale(0.0);
      hRK = RKPMBasisFunction(RKPMBasisOrder, xi, xi);
      phi = PHI(0.0);
      M.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, phi, hRK, hRK, 1.0);
      
      int mMatrixneighborhoodListIndex = neighborListIndex;
      const double* mbondDamage = bondDamage;
      for(int iNID=0; iNID<numNeighbors; ++iNID, ++mbondDamage)
      {
          neighborIndex = localNeighborList[mMatrixneighborhoodListIndex++];
          xn[0] = xOverlap[3*neighborIndex];
          xn[1] = xOverlap[3*neighborIndex+1];
          xn[2] = xOverlap[3*neighborIndex+2];
          distanceKernel = sqrt( (xi[0]-xn[0])*(xi[0]-xn[0]) + (xi[1]-xn[1])*(xi[1]-xn[1]) + (xi[2]-xn[2])*(xi[2]-xn[2]) );
          hRK = RKPMBasisFunction(RKPMBasisOrder, xi, xn);
          double zeta = std::abs(distanceKernel/RKPMSupport);
          phi = PHI(zeta)*(1-*mbondDamage);
          M.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, phi, hRK, hRK, 1.0);
      }

      //This loop is to calculate shape functions at the nodes in the
      //neighborhood. This will use the M matrix calculated from the loop above
      
      Teuchos::SerialDenseMatrix<int, double> hRK0(RKDimension, 1);
      hRK0.scale(0.0);
      hRK0(0,0) = 1.0;
      
      Teuchos::SerialQRDenseSolver<int, double> solver;
      solver.setMatrix(Teuchos::rcp(&M, false));
      solver.setVectors(Teuchos::rcp(&b, false), Teuchos::rcp(&hRK0, false));
      solver.solve();

      hRK = RKPMBasisFunction(RKPMBasisOrder, xi, xi);
      phi = PHI(0.0);
      Teuchos::SerialDenseMatrix<int, double> shapeFunction(1, 1);
      shapeFunction.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, phi, hRK, b, 0.0);

      OwnRKPMShapeFunctionValues[p] = shapeFunction(0,0);
      double consistency = 0.0;
      consistency += shapeFunction(0,0);
      for(int iNID=0 ; iNID<numNeighbors ; ++iNID, ++bondDamage)
      {
          neighborIndex = localNeighborList[neighborListIndex++];
          xn[0] = xOverlap[3*neighborIndex];
          xn[1] = xOverlap[3*neighborIndex+1];
          xn[2] = xOverlap[3*neighborIndex+2];
          distanceKernel = sqrt( (xi[0]-xn[0])*(xi[0]-xn[0]) + (xi[1]-xn[1])*(xi[1]-xn[1]) + (xi[2]-xn[2])*(xi[2]-xn[2]) );
          hRK = RKPMBasisFunction(RKPMBasisOrder, xi, xn);
          double zeta = std::abs(distanceKernel/RKPMSupport);
          phi = PHI(zeta)*(1-*bondDamage);
          shapeFunction.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, phi, hRK, b, 0.0);
          BondRKPMShapeFunctionValues[rkpmFunctionIndex++] = shapeFunction(0,0);
          consistency += shapeFunction(0,0);
          //std::cout << shapeFunction(0,0) << " ";
      }
      //std::cout << "|" << consistency << " ";
      //std::cout <<  std::endl;
      
      //Test if the RKPM ShapeFunction statisfies 0th order consistency 
      if (std::abs(consistency-1.0)>1e-10)
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Calculated RKPM Shape Functions do not satisfy 0th order consistency condition\n");
    }
    //std::cout << "computeRKPM" << std::endl;
    //*/
}

Teuchos::SerialDenseMatrix<int, double> RKPMBasisFunction(int RKPMBasisOrder, double *xi, double *xn) 
{

    //int RKDimension_total = RKPMBasisOrder+1; //1D
    //int RKDimension_total = (RKPMBasisOrder+1)(RKPMBasisOrder+2)/2; //2D
    int RKDimension_total = (RKPMBasisOrder+1)*(RKPMBasisOrder+2)*(RKPMBasisOrder+3)/6; //3D
    
    Teuchos::SerialDenseMatrix<int, double> RK(RKDimension_total, 1);
    int RKdimension = 0;
    int xOrder = RKPMBasisOrder;
    int yOrder = RKPMBasisOrder;
    int zOrder = RKPMBasisOrder;

    for(int kz=0; kz<=zOrder; ++kz)
    {
        for(int jy=0; jy<=yOrder; ++jy)
        {
            for(int ix=0; ix<=xOrder; ++ix)
            {
                if((ix+jy+kz) <= RKPMBasisOrder && RKdimension<RKDimension_total)
                {
                    RK(RKdimension, 0) = pow((xn[0]-xi[0]), ix) * pow((xn[1]-xi[1]), jy) * pow((xn[2]-xi[2]), kz);
                    RKdimension++;
                }
            }
        }
    }
    RK(0,0) = 1.0;
    return RK;
}
/*
double RKPMCubicSplineKernel(double a, double distanceKernel) 
{
    double phi;

    double z = std::abs(distanceKernel)/a;
    if (z >= 0.0 & z <= 0.5)
        phi = 2.0/3.0 - 4.0*pow(z, 2) + 4.0*pow(z, 3);
    else if (z > 0.5 && z <= 1.0)
        phi = 4.0/3.0 * pow((1.0- z), 3);
    else
        phi = 0.0;
    return phi;
}
*/
void applyShapeFunction
( 
        const double* xOverlap, 
        double* yOverlap, 
        int myNumPoints, 
        const int* localNeighborhoodList, 
        const double* OwnRKPMShapeFunctionValues, 
        const double* BondRKPMShapeFunctionValues
) 
{
    ///*
    int numNeighbors, neighborId, neighborhoodListIndex(0), rkpmShapeFunctionIndex(0);
    for(int iID=0; iID<myNumPoints; ++iID)
    {
        double consistency = 0.0;
        double RKPMShapeFunction = 0.0;
        RKPMShapeFunction = OwnRKPMShapeFunctionValues[iID];
        consistency += RKPMShapeFunction;
        yOverlap[iID*3+0] = xOverlap[iID*3+0] + (yOverlap[iID*3+0] - xOverlap[iID*3+0]) * RKPMShapeFunction;
        yOverlap[iID*3+1] = xOverlap[iID*3+1] + (yOverlap[iID*3+1] - xOverlap[iID*3+1]) * RKPMShapeFunction;
        yOverlap[iID*3+2] = xOverlap[iID*3+2] + (yOverlap[iID*3+2] - xOverlap[iID*3+2]) * RKPMShapeFunction;
        numNeighbors = localNeighborhoodList[neighborhoodListIndex++];
        for(int iNID=0 ; iNID<numNeighbors ; ++iNID)
        {
            neighborId = localNeighborhoodList[neighborhoodListIndex++];
            RKPMShapeFunction = BondRKPMShapeFunctionValues[rkpmShapeFunctionIndex++];
            consistency += RKPMShapeFunction;
            yOverlap[iID*3+0] += (yOverlap[neighborId*3+0] - xOverlap[neighborId*3+0]) * RKPMShapeFunction; 
            yOverlap[iID*3+1] += (yOverlap[neighborId*3+1] - xOverlap[neighborId*3+1]) * RKPMShapeFunction; 
            yOverlap[iID*3+2] += (yOverlap[neighborId*3+2] - xOverlap[neighborId*3+2]) * RKPMShapeFunction; 
        } 
        //cout << consistency << " ";
        
        //Test if the RKPM ShapeFunction statisfies 0th order consistency 
        if (std::abs(consistency-1.0)>1e-10)
        {
          std::cout << consistency << " ";
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** RKPM Shape Functions do not satisfy 0th order consistency condition\n");
        }

   } //*/
   //std::cout << "applyRKPM" << std::endl;

}

}
