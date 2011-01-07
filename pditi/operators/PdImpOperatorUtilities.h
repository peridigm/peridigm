/*
 * PimpUtilitities.h
 *
 *  Created on: Apr 3, 2010
 *      Author: awesome
 */

#ifndef PIMP_OPERATOR_UTILITIES_H_
#define PIMP_OPERATOR_UTILITIES_H_
#include <math.h>
#include "Field.h"

namespace PdImp {

/*
 * INCOMING Arrays must have a minimum dimension of 3
 * bondIP = xP - xI
 */
inline void BOND(const double *xI, const double* xP, double *bondIP){
	* bondIP =    * xP    - * xI;
	*(bondIP+1) = *(xP+1) - *(xI+1);
	*(bondIP+2) = *(xP+2) - *(xI+2);
}

inline void UPDATE_GEOMETRY(const double *x, const double* u, double *y){
	* y    = * x     + * u    ;
	*(y+1) = *(x+1)  + *(u+1) ;
	*(y+2) = *(x+2)  + *(u+2) ;
}

inline void UPDATE_GEOMETRY(const double *x, const double* u, double *y, int numPoints){
	for(int p=0;p<numPoints;p++,x+=3,u+=3,y+=3){
		* y    = * x     + * u    ;
		*(y+1) = *(x+1)  + *(u+1) ;
		*(y+2) = *(x+2)  + *(u+2) ;
	}
}

inline double MAGNITUDE(const double *x){
	double x1 = * x;
	double x2 = *(x+1);
	double x3 = *(x+2);
	return sqrt(x1*x1+x2*x2+x3*x3);
}

/**
 *@param x -- this vector is normalized
 *@return magnitude of x
 */
inline double NORMALIZE(double *x){
	double norm = MAGNITUDE(x);
	* x    = *(x  )/norm;
	*(x+1) = *(x+1)/norm;
	*(x+2) = *(x+2)/norm;
	return norm;
}

inline double DOT(const double *u, const double *v){
	double dot = 0;
	dot += (*(u))   * (*(v));
	dot += (*(u+1)) * (*(v+1));
	dot += (*(u+2)) * (*(v+2));
	return dot;
}

/**
 * Forms matrix k from the tensor product of x with y in that order
 * @param x -- vector of length 3
 * @param y -- vector of length 3
 * @param k -- vector of length 9
 * @return k -- column major storage
 */
void TENSOR_PRODUCT(const double *x, const double *y, double *k);

inline void SET(double* const start, const double* const end, double value){
	for(double *i = start; i!=end; i++)
			*i = value;
}

inline double* SUMINTO(const double* const start, const double* const end, double* intoMe){
	double *out = intoMe;
	for(const double *i = start; i!=end; i++, out++)
			*out += *i;
	return intoMe;
}

inline double* SUBTRACTINTO(const double* const start, const double* const end, double* intoMe){
	double *out = intoMe;
	for(const double *i = start; i!=end; i++, out++)
			*out -= *i;
	return intoMe;
}

inline double* COPY(const double* const start, const double* const end, double* intoMe){
	double *out = intoMe;
	for(const double *i = start; i!=end; i++, out++)
			*out = *i;
	return intoMe;
}

inline double* SCALE_BY_VALUE(double* scaleMe, const double* const scaleMeEND, double value){
	for(double *x = scaleMe; x!=scaleMeEND; x++)
		*x *= value;
	return scaleMe;
}

void PRINT_3x3MATRIX(const double * const k, std::ostream& out);

/**
 * @param owned point coordinates
 * @return displacement at points that create pure shear
 */
Field_NS::Field<double> getPureShearXY(const Field_NS::Field<double>& X);

} // Namespace PdImpOperators

#endif /* PIMP_OPERATOR_UTILITIES_H_ */
