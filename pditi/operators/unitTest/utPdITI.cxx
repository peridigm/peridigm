/*
 * utPdITI.cxx
 *
 *  Created on: Feb 22, 2011
 *      Author: jamitch
 */

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>
#include <map>
#include <vector>
#include <set>

#include "Array.h"
#include "utPdITI.h"
#include "../DirichletBcSpec.h"


#include "Epetra_BlockMap.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_RowMatrix.h"
#include <Epetra_FEVbrMatrix.h>

using std::vector;
using std::set;
using PdImp::DirichletBcSpec;
using UTILITIES::Array;
using namespace PdImp;

const int vectorNDF=3;


namespace utPdITI {

shared_ptr<Epetra_CrsGraph> getGraph(shared_ptr<RowStiffnessOperator>& jacobian){
	const Epetra_BlockMap& rowMap   = jacobian->getRowMap();
	const Epetra_BlockMap& colMap = jacobian->getColMap();

	/*
	 * Epetra Graph
	 */
	Array<int> numCols = jacobian->getNumColumnsPerRow();
	Epetra_CrsGraph *graph = new Epetra_CrsGraph(Copy,rowMap,colMap,numCols.get());
	for(int row=0;row<rowMap.NumMyElements();row++){
		Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numCol = rowLIDs.get_size();
		int *cols = rowLIDs.get();
		if(0!=graph->InsertMyIndices(row,numCol,cols)){
			std::string message("utPdITI::graph->InsertMyIndices(row,numCol,cols)\n");
			message += "\t0!=graph->FillComplete()";
			throw std::runtime_error(message);
		}
	}

	if(0!=graph->FillComplete()){
		std::string message("utPdITI::getGraph(shared_ptr<RowStiffnessOperator>& jacobian)\n");
		message += "\t0!=graph->FillComplete()";
		throw std::runtime_error(message);
	}
	return shared_ptr<Epetra_CrsGraph>(graph);
}

IsotropicHookeSpec getMaterialSpec(double e, double nu) {
	YoungsModulus youngsModulus = IsotropicHookeSpec::youngsModulus(e);
	PoissonsRatio poissonsRatio = IsotropicHookeSpec::poissonsRatio(nu);
	return IsotropicHookeSpec(youngsModulus,poissonsRatio);
}

shared_ptr<Epetra_RowMatrix>
getOperator
(
		const vector<shared_ptr<StageComponentDirichletBc> >& bcArray,
		shared_ptr<Epetra_CrsGraph>& graphPtr,
		shared_ptr<RowStiffnessOperator>& jacobian
)
{
	std::cout << "Begin jacobian calculation\n";
	const Epetra_BlockMap& rowMap   = jacobian->getRowMap();

	/*
	 * Loop over Bc's and create set of ids that can be searched
	 */

	vector<std::set<int> > bcPointIds(bcArray.size());
	{
		vector<shared_ptr<StageComponentDirichletBc> >::const_iterator bcIter = bcArray.begin();
		for(int i=0;bcIter != bcArray.end(); i++,bcIter++){
			StageComponentDirichletBc* stageComponentPtr = bcIter->get();
			const DirichletBcSpec& spec = stageComponentPtr->getSpec();
			const Array<int>& ids = spec.getPointIds();
			bcPointIds[i]= std::set<int>(ids.get(),ids.get()+ids.get_size());
		}
	}

	/*
	 * Create a searchable set for bc
	 */

	/*
	 * Epetra Matrix
	 * PERHAPS the 'operator' can keep its own copy of Graph, then this
	 * can be a 'View'
	 */

	Epetra_FEVbrMatrix *m = new Epetra_FEVbrMatrix(Copy,*(graphPtr.get()));

	Epetra_SerialDenseMatrix k;
	k.Shape(vectorNDF,vectorNDF);
	for(int row=0;row<rowMap.NumMyElements();row++){
		Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numCol = rowLIDs.get_size();
		int *cols = rowLIDs.get();

		if(0!=m->BeginReplaceMyValues(row,numCol,cols)){
			std::string message("utPdITI::getOperator(bcArray,graphPtr,jacobian)\n");
			message += "\t0!=m->BeginReplaceMyValues(row,numCol,cols)";
			throw std::runtime_error(message);
		}

		/*
		 * loop over columns in row and submit block entry
		 */
		Array<double> actualK = jacobian->computeRowStiffness(row, rowLIDs);


		/*
		 * 1) Zero out row as necessary due to BC's
		 * 2) Assemble into matrix -- zero columns as necessary
		 * 3)
		 */
		vector<std::set<int> >::iterator pointSetIter = bcPointIds.begin();
		vector<shared_ptr<StageComponentDirichletBc> >::const_iterator bcIter = bcArray.begin();
		for(;bcIter != bcArray.end(); bcIter++, pointSetIter++){
			const std::set<int>::const_iterator bcIdsEnd = pointSetIter->end();

			/*
			 * Get components to be applied
			 */
			StageComponentDirichletBc* stageComponentPtr = bcIter->get();
			const DirichletBcSpec& spec = stageComponentPtr->getSpec();
			vector<DirichletBcSpec::ComponentLabel> components = spec.getComponents();

			/*
			 * Search for row in bcIds
			 */
			// if this is true, then this row is a bc row
			if(bcIdsEnd != pointSetIter->find(row)) {

				/*
				 * 1) This row has a bc; need to zero out in stiffness
				 * 2) Watch for diagonal: place "1" on the diagonal
				 */

				/*
				 * Create array of pointers to each row/component that BC is applied
				 */
				vector<double*> rowPtrs(components.size(), NULL);
				vector<double*> diagPtrs(components.size(), NULL);
				for(std::size_t r=0;r<components.size();r++)
					rowPtrs[r]= actualK.get()+components[r];

				for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++){
					int col = *colPtr;

					for(std::size_t r=0;r<components.size();r++){
						/*
						 * 0) Save diagonal location
						 * 1) Set row element to zero
						 * 2) Move to next column
						 * 3) Note that there are 3 columns in 3x3 matrix
						 */

						/*
						 * Save location of diagonal for this component
						 */
						diagPtrs[r] = rowPtrs[r] + 3 * components[r];

						/*
						 * Zero out row corresponding with component
						 * Increment row pointer to next column
						 */
						*(rowPtrs[r])=0; rowPtrs[r]+=3;
						*(rowPtrs[r])=0; rowPtrs[r]+=3;
						*(rowPtrs[r])=0; rowPtrs[r]+=3;

					}
					/*
					 * If this is true, then we are on the diagonal
					 * Need to put "1" on the strict diagonal by component
					 * Put -1 because later the whole matrix is negated
					 */
					if(row==col) {
						for(std::size_t r=0;r<components.size();r++){
							*(diagPtrs[r]) = -1.0;
						}

					}

				}

			} else {
				/*
				 * This is not a bc row but we still have to check for columns that may have a bc assigned to them
				 * In this case, we just have to zero the column
				 */
				double *kPtr=actualK.get();
				for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++,kPtr+=9){
					int col = *colPtr;
					if(bcIdsEnd == pointSetIter->find(col)) continue;
					/*
					 * We have a column that must be zero'd
					 */
					/*
					 * Zero out column for each component that is applied
					 */
					double *stiffPtr=0;
					for(std::size_t r=0;r<components.size();r++){
						stiffPtr = kPtr+3*components[r];
						for(int r=0;r<3;r++)
							stiffPtr[r]=0;
					}

				}

			}

		}

		/*
		 * Now just populate the matrix
		 */
		double *kPtr = actualK.get();
		for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++){

			/*
			 * Fill 'k'
			 */
			for(int c=0;c<3;c++){
				double *colK = k[c];
				for(int r=0;r<3;r++,kPtr++)
					colK[r] = - *kPtr;
			}

			if(0!=m->SubmitBlockEntry(k)){
				std::string message("utPdITI::getOperator(bcArray,graphPtr,jacobian)\n");
				message += "\t0!=m->SubmitBlockEntry(k)";
				throw std::runtime_error(message);
			}

		}


		/*
		 * Finalize this row
		 */
		if(0!=m->EndSubmitEntries()){
			std::string message("utPdITI::getOperator(bcArray,graphPtr,jacobian)\n");
			message += "\t0!=m->EndSubmitEntries()";
			throw std::runtime_error(message);
		}

	}

	if(0!=m->FillComplete()){
		std::string message("utPdITI::getOperator(bcArray,graphPtr,jacobian)\n");
		message += "\t0!=m->FillComplete()";
		throw std::runtime_error(message);
	}


	std::cout << "END jacobian calculation\n";
	return shared_ptr<Epetra_RowMatrix>(m);
}

}
