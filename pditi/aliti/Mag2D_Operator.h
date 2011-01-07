/*
 * Mag2D_Operator.h
 *
 *  Created on: Nov 18, 2010
 *      Author: jamitch
 */

#ifndef MAG2D_OPERATOR_H_
#define MAG2D_OPERATOR_H_

#include <iostream>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface_Required.H" // base class
#include "NOX_Epetra_Interface_Jacobian.H" // base class


class Mag2D_Operator :
public NOX::Epetra::Interface::Required,
public NOX::Epetra::Interface::Jacobian
{
public:
	Mag2D_Operator();
	~Mag2D_Operator();

	//! Compute and return F
	bool computeF(const Epetra_Vector& x, Epetra_Vector& FVec,FillType flag = Residual);

	//! Compute an explicit Jacobian
	bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

};



#endif /* MAG2D_OPERATOR_H_ */
