/*
 * PimpOperatorUtilities.cxx
 *
 *  Created on: Apr 13, 2010
 *      Author: jamitch
 */
#include "PdITI_Utilities.h"
#include <ostream>

namespace PdITI {

using Field_NS::Field;

void PRINT_3x3MATRIX(const double * const k, std::ostream& out){
	const int lda=3;
	for(int r=0;r<lda;r++){
		out <<  std::scientific << "\t" << *(k+r) << " " << *(k+lda+r) << " " << *(k+2*lda+r) << std::endl;
	}
}

void TENSOR_PRODUCT(const double *x, const double *y, double *k){
	double x1 = * x;
	double x2 = *(x+1);
	double x3 = *(x+2);
	double y1 = * y;
	double y2 = *(y+1);
	double y3 = *(y+2);

	int lda = 3;
	// Column 1
	* k    = x1 * y1;
	*(k+1) = x2 * y1;
	*(k+2) = x3 * y1;
	// Column 2
	k += lda;
	* k    = x1 * y2;
	*(k+1) = x2 * y2;
	*(k+2) = x3 * y2;
	// Column 3
	k += lda;
	* k    = x1 * y3;
	*(k+1) = x2 * y3;
	*(k+2) = x3 * y3;
}

Field<double> getPureShearXY(const Field<double>& X){
	std::size_t numPoints = X.get_num_points();
	Field<double> uOwnedField = Field<double>(Field_NS::DISPL3D,numPoints);
	double *u = uOwnedField.get();
	const double *x = X.get();

	double gamma=2.0;
	for(std::size_t i=0;i<numPoints;i++){
		int p=3*i;
		u[p]=gamma*x[p+1];
		u[p+1]=0;
		u[p+2]=0;
	}
	return uOwnedField;
}

}
