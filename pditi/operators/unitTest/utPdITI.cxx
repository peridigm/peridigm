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


#include "Epetra_BlockMap.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_RowMatrix.h"
#include <Epetra_FEVbrMatrix.h>

using std::vector;
using std::set;
using PdImp::DirichletBcSpec;
using Field_NS::Field;
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
        const Field<char> bcMaskFieldOverlap,
        shared_ptr<RowStiffnessOperator>& jacobian
)
{
    std::cout << "Begin jacobian calculation utPdITI::getOperator(...)\n";
    const Epetra_BlockMap& rowMap   = jacobian->getRowMap();
    const Epetra_BlockMap& colMap   = jacobian->getColMap();
    Array<int> numColPerRow;
    numColPerRow.deep_copy(jacobian->getNumColumnsPerRow());
    Epetra_FEVbrMatrix *m = new Epetra_FEVbrMatrix(Copy,rowMap,numColPerRow.get());

    const char *bcMask = bcMaskFieldOverlap.get();
    Epetra_SerialDenseMatrix k;
    k.Shape(vectorNDF,vectorNDF);
    /*
     * DEBUG PRINTING
     */
//    std::stringstream streamRows,streamNumRows,streamNumCols, streamBcMask;
//    streamNumRows << "int numRows(" << rowMap.NumMyElements() << ");";
//    streamNumCols << "int numCols[]="<< "{";
//
//    for(std::size_t row=0;row<rowMap.NumMyElements();row++){
//    	if(0==row){
//    		streamRows << "int* rowPtrPID" << rowMap.Comm().MyPID() << "[]={";
//    		streamBcMask << "char bcMaskRowPID" << rowMap.Comm().MyPID() << "[]={";
//    	}
//    	if(row%10==0){
//    		streamRows << "\n";
//    		streamNumCols << "\n";
//    		streamBcMask << "\n";
//    	}
//    	streamRows << "rowGID" << rowMap.GID(row);
//    	int b = bcMask[row];
//    	streamBcMask << b;
//    	Array<int> rowLIDs = jacobian->getColumnLIDs(row);
//    	streamNumCols << rowLIDs.get_size();
//    	if(row!=rowMap.NumMyElements()-1){
//    		streamNumCols << ", ";
//    		streamRows << ", ";
//    		streamBcMask << ", ";
//    	} else
//    	{
//    		streamNumCols << "};\n";
//    		streamRows << "};\n";
//    		streamBcMask << "};\n";
//    	}
//    }
//    std::string strNumRows=streamNumRows.str();
//    std::string strNumCols=streamNumCols.str();
//    std::cout << strNumRows << std::endl;
//    std::cout << streamRows.str() << std::endl;
//    std::cout << strNumCols << std::endl;
//    std::cout << streamBcMask.str() << std::endl;
    /*
     * DEBUG PRINTING
     */


    for(int row=0;row<rowMap.NumMyElements();row++){
        int rowGID = rowMap.GID(row);
        Array<int> rowLIDs = jacobian->getColumnLIDs(row);
        std::size_t numCol = rowLIDs.get_size();

        /*
         * Hack to convert rowLIDs over to rowGIDs; If this approach
         * to assembly works, then eliminate the rowLIDs call and create an
         * analogous call 'getColumnGIDs(lowRow)
         */
        Array<int> colGIDs(numCol);
        for(size_t c=0;c<numCol;c++)
            colGIDs[c]=colMap.GID(rowLIDs[c]);

        if(0!=m->BeginInsertGlobalValues(rowGID,numCol,colGIDs.get())){
            std::string message("utPdITI::getOperator_NEWER(bcArray,jacobian)\n");
            message += "\t0!=m->BeginReplaceGlobalValues(rowGID,numCol,colGIDs.get())";
            throw std::runtime_error(message);
        }

    	/*
    	 * Print rows
    	 */

//        std::stringstream streamRow, streamColBcMask;
//        streamRow << "int rowGID" << rowGID << "[]={";
//        streamColBcMask << "int bcColMaskGID" << rowGID << "[]={";
//        for(int i=0;i<colGIDs.get_size();i++){
//        	if(i%10==0)
//        		streamRow << "\n";
//        	streamRow << *(colGIDs.get()+i);
//        	if(i!=colGIDs.get_size()-1)
//        		streamRow << ", ";
//        	if(i==colGIDs.get_size()-1)
//        		streamRow << "};\n";
//        }
//        std::cout << streamRow.str() << std::endl;

    	/*
    	 * End print rows
    	 */


        /*
          * Dirichlet BC for this row
          */
 		char rowMask = bcMask[row];


		/*
		 * Create array of pointers to each row/component that BC is applied
		 */
		vector<DirichletBcSpec::ComponentLabel> rowComponentBcs = getComponents(rowMask);
		vector<double*> rowPtrs(rowComponentBcs.size(), NULL);
		vector<double*> diagPtrs(rowComponentBcs.size(), NULL);
		Array<double> actualK = jacobian->computeRowStiffness(row, rowLIDs);

 		/*
 		 * Submit block row entries based upon GID
 		 * GIDs for each column in row
 		 */
        int *cols = colGIDs.get();
		double *colRoot = actualK.get();
		size_t dum=0;
		for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++, colRoot+=9, dum++){
			int col = *colPtr;

			/*
			 * Handle BCs for row
			 */
			for(std::size_t r=0;r<rowComponentBcs.size();r++){
				rowPtrs[r] = colRoot + rowComponentBcs[r];
				/*
				 * 0) Save diagonal location
				 * 1) Set row element to zero
				 * 2) Move to next column
				 * 3) Note that there are 3 columns in 3x3 matrix :) DUH!
				 */

				/*
				 * Save location of diagonal for this component
				 */
				diagPtrs[r] = rowPtrs[r] + 3 * rowComponentBcs[r];

				/*
				 * Zero out row corresponding with component
				 * Increment row pointer to next column
				 */
				*(rowPtrs[r])=0; rowPtrs[r]+=3;
				*(rowPtrs[r])=0; rowPtrs[r]+=3;
				*(rowPtrs[r])=0; rowPtrs[r]+=3;

			}

			/*
			 * Now check column for dirichlet bc; here we will just zero everything out;
			 * If this column happens to be the diagonal, then we will fix later;
			 */
			char colMask = bcMask[colMap.LID(col)];
			/*
			 * DEBUG PRINTING
			 */
//			if(dum%10==0){
//				streamColBcMask << "\n";
//			}
//			{
//				int b = colMask;
//				streamColBcMask << b;
//			}
//			if(dum!=numCol-1)
//				streamColBcMask << ", ";
//			if(dum==numCol-1)
//				streamColBcMask << "};\n";
			/*
			 * END DEBUG PRINTING
			 */
			vector<DirichletBcSpec::ComponentLabel> colComponentBcs = getComponents(colMask);
			double *stiffPtr(0);
			for(std::size_t c=0;c<colComponentBcs.size();c++){
				stiffPtr = colRoot + 3 * colComponentBcs[c];
				for(int r=0;r<3;r++)
					stiffPtr[r]=0;
			}


			/*
			 * If this is true, then we are on the diagonal
			 * Need to put "1" on the strict diagonal by component
			 * Put -1 because later the whole matrix is negated
			 */
			if(rowGID==col) {
				for(std::size_t r=0;r<rowComponentBcs.size();r++){
					*(diagPtrs[r]) = -1.0;
				}
			}


		}
		/*
		 * DEBUG PRINTING
		 */
//		std::cout << streamColBcMask.str() << std::endl;
		/*
		 * END DEBUG PRINTING
		 */

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
				std::string message("utPdITI::getOperator(bcArray,jacobian)\n");
				message += "\t0!=m->SubmitBlockEntry(k)";
				throw std::runtime_error(message);
			}

		}

		/*
		 * Finalize this row
		 */
		if(0!=m->EndSubmitEntries()){
			std::string message("utPdITI::getOperator(bcArray,jacobian)\n");
			message += "\t0!=m->EndSubmitEntries()";
			throw std::runtime_error(message);
		}

	}

    int err = m->FillComplete(rowMap,rowMap);
	if(0!=err){
		std::string message("utPdITI::getOperator(bcArray,jacobian)\n");
		message += "\t0!=m->FillComplete()";
		throw std::runtime_error(message);
	}

	std::cout << "END jacobian calculation utPdITI::getOperator(...)\n";
    return shared_ptr<Epetra_RowMatrix>(m);

}

vector<DirichletBcSpec::ComponentLabel> getComponents(char mask) {

	vector<DirichletBcSpec::ComponentLabel> c;
	if(1 & mask)
		c.push_back(DirichletBcSpec::X);
	if(2 & mask)
		c.push_back(DirichletBcSpec::Y);
	if(4 & mask)
		c.push_back(DirichletBcSpec::Z);

	return c;
}


}
