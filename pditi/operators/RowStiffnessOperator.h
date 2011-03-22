/*
 * RowStiffnessOperator.h
 *
 *  Created on: Mar 30, 2010
 *      Author: jamitch
 */

#ifndef ROWSTIFFNESSOPERATOR_H_
#define ROWSTIFFNESSOPERATOR_H_
#include "Array.h"

class Epetra_BlockMap;

class RowStiffnessOperator {
public:
	virtual ~RowStiffnessOperator() {};
	virtual int getNumRows() const = 0;
	virtual UTILITIES::Array<int> getColumnLIDs(int localRowID) const = 0;
	virtual const UTILITIES::Array<int>& getNumColumnsPerRow() const = 0;
	virtual const UTILITIES::Array<double>& computeRowStiffness(int localRowID, UTILITIES::Array<int> rowGIDs) = 0;
	virtual const Epetra_BlockMap& getRowMap() const = 0;
	virtual const Epetra_BlockMap& getColMap() const = 0;
};


#endif /* ROWSTIFFNESSOPERATOR_H_ */
