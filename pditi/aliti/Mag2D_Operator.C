/*
 * Mag2D_Operator.C
 *
 *  Created on: Nov 18, 2010
 *      Author: jamitch
 */

#include "Mag2D_Operator.h"

Mag2D_Operator::Mag2D_Operator()  {}
Mag2D_Operator::~Mag2D_Operator() {}

bool Mag2D_Operator::computeF(const Epetra_Vector& x, Epetra_Vector& FVec,FillType flag) {
	return false;
}

bool Mag2D_Operator::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) {
	return false;
}
