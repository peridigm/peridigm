/*
 * utCrackOpeningDemo.cxx
 *
 *  Created on: Feb 11, 2011
 *      Author: jamitch
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_ALTERNATIVE_INIT_API
//#include <boost/test/unit_test.hpp>
//#include <boost/test/parameterized_test.hpp>
#include "Sortable.h"
#include "Array.h"
#include "quick_grid/QuickGrid.h"
#include "NeighborhoodList.h"
#include "OverlapDistributor.h"
#include "BondFilter.h"
#include "PdZoltan.h"
#include "vtk/Field.h"
#include "vtk/PdVTK.h"
#include "../../PdImpMaterials.h"
#include "../../PdITI_Operator.h"
#include "../../PdITI_Utilities.h"
#include "../../DirichletBcSpec.h"
#include "../../BodyLoadSpec.h"
#include "../../StageFunction.h"
#include "../../Loader.h"
#include "../../StageComponentDirichletBc.h"
#include "../../ComponentDirichletBcSpec.h"
#include "../../IsotropicElasticConstitutiveModel.h"
#include "../../ConstitutiveModel.h"
#include "../utPdITI.h"
#include "PdutMpiFixture.h"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <set>
#include <Epetra_FEVbrMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <utility>
#include <map>
#include <vector>
#include <cstdlib>
#include <time.h>


using UTILITIES::CartesianComponent;
using UTILITIES::Array;
using PdITI::ConstitutiveModel;
using PdITI::IsotropicElasticConstitutiveModel;
using PdImp::StageComponentDirichletBc;
using PdImp::StageDirichletBc;
using PdImp::DirichletBcSpec;
using PdImp::ComponentDirichletBcSpec;
using PdImp::StageFunction;
using namespace PdBondFilter;
using namespace Pdut;
using namespace Field_NS;
using std::tr1::shared_ptr;

const int vectorNDF=3;
static size_t myRank;
static size_t numProcs;

/*
 * This should be even so that the crack plane lies between to rows of points
 */
const size_t nx = 20;
const size_t ny = 20;
const double xStart = -2.5;
const double xLength = 5.0;
const double yStart = -2.5;
const double yLength = 5.0;
const double zStart = -0.5;
const double zLength = 1.0;
static double xMax=xStart+xLength;
static double xMin=xStart;
const size_t nz = nx* (zLength / xLength);
const size_t numCells = nx*ny*nz;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const double dx = xSpec.getCellSize();
const double dy = ySpec.getCellSize();
const double dz = zSpec.getCellSize();
const double _cellVolume = dx*dy*dz;

/*
 * Young's Modulus (MPa)
 */
const double E = 68.9e3;

/*
 * Poisson's ratio
 */
const double nu = 0.0;

/*
 * Density of aluminum g/mm^3
 */
const double rho = 2.7e-3;

/*
 * Horizon
 */
const double horizon=1.1*sqrt( (3.0*dx)*(3.0*dx) );

/*
 * Function prototypes in this file
 */
FinitePlane getYZ_CrackPlane();
vector<DirichletBcSpec::ComponentLabel> getComponents(char mask);


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
			std::string message("graph->InsertMyIndices(row,numCol,cols)\n");
			message += "\t0!=graph->FillComplete()";
			throw std::runtime_error(message);
		}
	}

	if(0!=graph->FillComplete()){
		std::string message("getGraph(shared_ptr<RowStiffnessOperator>& jacobian)\n");
		message += "\t0!=graph->FillComplete()";
		throw std::runtime_error(message);
	}
	return shared_ptr<Epetra_CrsGraph>(graph);
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
			bcPointIds[i]= std::set<int>(ids.get(),ids.end());
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

shared_ptr<Epetra_RowMatrix>
getOperator_NEW
(
		const Field<char> bcMaskFieldOverlap,
		shared_ptr<Epetra_CrsGraph>& graphPtr,
		shared_ptr<RowStiffnessOperator>& jacobian
)
{
	std::cout << "Begin jacobian calculation\n";
	const Epetra_BlockMap& rowMap   = jacobian->getRowMap();

	/*
	 * Epetra Matrix
	 * PERHAPS the 'operator' can keep its own copy of Graph, then this
	 * can be a 'View'
	 */

	Epetra_FEVbrMatrix *m = new Epetra_FEVbrMatrix(Copy,*(graphPtr.get()));

	const char *bcMask = bcMaskFieldOverlap.get();
	Epetra_SerialDenseMatrix k;
	k.Shape(vectorNDF,vectorNDF);
	for(int row=0;row<rowMap.NumMyElements();row++){
		Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numCol = rowLIDs.get_size();

		std::vector< std::pair<int,int> > pairs(numCol);
		{
			int *cols = rowLIDs.get();
			for(std::size_t i=0;i<pairs.size();i++,cols++){
				pairs[i] = std::make_pair(*cols,i);
			}
		}
		std::map<int,int> map(pairs.begin(),pairs.end());


		int *cols = rowLIDs.get();
		char rowMask = bcMask[row];


		if(0!=m->BeginReplaceMyValues(row,numCol,cols)){
			std::string message("utPdITI::getOperator(bcArray,graphPtr,jacobian)\n");
			message += "\t0!=m->BeginReplaceMyValues(row,numCol,cols)";
			throw std::runtime_error(message);
		}

		/*
		 * Create array of pointers to each row/component that BC is applied
		 */
		vector<DirichletBcSpec::ComponentLabel> rowComponentBcs = getComponents(rowMask);
		vector<double*> rowPtrs(rowComponentBcs.size(), NULL);
		vector<double*> diagPtrs(rowComponentBcs.size(), NULL);
		Array<double> actualK = jacobian->computeRowStiffness(row, rowLIDs);

		/*
		 * Zero out row as necessary due to BC's
		 */
		double *colRoot;
		for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++){
			int col = *colPtr;
			char colMask = bcMask[col];

			/*
			 * Handle BCs for row
			 */
			colRoot = actualK.get() + 9 * map[col];
			for(std::size_t r=0;r<rowComponentBcs.size();r++){
				rowPtrs[r]= colRoot + rowComponentBcs[r];
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
			if(row==col) {
				for(std::size_t r=0;r<rowComponentBcs.size();r++){
					*(diagPtrs[r]) = -1.0;
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

shared_ptr<Epetra_RowMatrix>
getOperator_NEWER
(
        const Field<char> bcMaskFieldOverlap,
        shared_ptr<Epetra_CrsGraph>& graphPtr,
        shared_ptr<RowStiffnessOperator>& jacobian
)
{
    std::cout << "Begin jacobian calculation NEWER\n";
    const Epetra_BlockMap& rowMap   = jacobian->getRowMap();
    const Epetra_BlockMap& colMap   = jacobian->getColMap();
    Array<int> numColPerRow;
    numColPerRow.deep_copy(jacobian->getNumColumnsPerRow());
    Epetra_FEVbrMatrix *m = new Epetra_FEVbrMatrix(Copy,rowMap,numColPerRow.get());

    const char *bcMask = bcMaskFieldOverlap.get();
    Epetra_SerialDenseMatrix k;
    k.Shape(vectorNDF,vectorNDF);
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
            std::string message("utPdITI::getOperator_NEWER(bcArray,graphPtr,jacobian)\n");
            message += "\t0!=m->BeginReplaceGlobalValues(rowGID,numCol,colGIDs.get())";
            throw std::runtime_error(message);
        }

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
		for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++, colRoot+=9){
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

/*
 * This demonstrates how the first and second coordinate
 * along an axis are computed
 */
//const double x0 = xStart+xSpec.getCellSize()/2.0;
//const double x1 = x0 + xSpec.getCellSize();
//const double y0 = yStart+ySpec.getCellSize()/2.0;
//const double y1 = y0 + ySpec.getCellSize();
//const double z0 = zStart+zSpec.getCellSize()/2.0;
//const double z1 = z0 + zSpec.getCellSize();


QUICKGRID::QuickGridData getGrid() {

	if(0==myRank){
		cout << "Creating and load balancing mesh..." << endl;
	}

	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData gridData =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);
	gridData=PDNEIGH::getLoadBalancedDiscretization(gridData);

	/*
	 * Lower left hand corner of crack plane when viewing down
	 * normal in the +dir
	 */
	const double x0 = xStart+xLength/2;
	const double y0 = yStart;
	const double z0 = zStart;

	if(0==myRank){
		cout << "\t\tDONE." << endl;
		cout << "Total number of points in mesh = " << gridData.globalNumPoints << endl;
		cout << "nx,ny,nz = " << nx << ", " << ny << ", "<< nz << ", "<< endl;
		cout << "x0,y0,z0 = " << x0 << ", " << y0 << ", "<< z0 << ", "<< endl;
	}

	return gridData;
}



FinitePlane getYZ_CrackPlane() {

	/*
	 * Crack runs along y-axis
	 * Crack length along y-axis is yLength/2
	 * Crack runs from bottom to top of plate in z-dir; a=zLength
	 */

	/*
	 * Lower left hand corner of crack plane when viewing down
	 * normal in the +dir
	 */
	const double x0 = xStart+xLength/2;
	const double y0 = yStart;
	const double z0 = zStart;

	/*
	 * normal is along -x-dir
	 */
	double n[3]; n[0]=-1.0;n[1]=0.0;n[2]=0.0;
	/*
	 * lower left corner of plane
	 */
	double r0[3]; r0[0]=x0; r0[1]=y0; r0[2]=z0;
	/*
	 * vector along bottom edge is in the +y-dir
	 */
	double ub[3]; ub[0]=0; ub[1]=1.0;ub[2]=0.0;
	/*
	 * b is length of crack along bottom edge
	 * a is length of crack along z-dir
	 */
	double b=yLength/2.0, a=zLength;
	return FinitePlane(n,r0,ub,b,a);
}

void crackOpeningDemo(){
	/*
	 * Get mesh and decomposition
	 */
	QUICKGRID::QuickGridData gridData = getGrid();
	/*
	 * Communicator
	 */
	shared_ptr<Epetra_Comm> comm = shared_ptr<Epetra_Comm>(new Epetra_MpiComm(MPI_COMM_WORLD));

	FinitePlane crackPlane=getYZ_CrackPlane();
//	shared_ptr<BondFilter> filterPtr = shared_ptr<BondFilter>(new BondFilterDefault());
	shared_ptr<BondFilter> filterPtr = shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	/*
	 * Material Properties
	 */
	IsotropicHookeSpec isotropicSpec = utPdITI::getMaterialSpec(E,nu);

	/*
	 * Create PdITI Operator
	 */
	PdITI::PdITI_Operator op(comm,list,gridData.cellVolume);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	op.addConstitutiveModel(fIntOperator);
	PDNEIGH::NeighborhoodList row_matrix_list = op.get_row_matrix_neighborhood();

	/*
	 * Get points for bc's
	 */
	/*
	 * Note that we are looking for a discrete number of points at end;
	 * Set the scale factor to just larger than an integer where
	 * the integer corresponds with the number of points to be included in
	 * the boundary conditions
	 */
	double scaleFactor=2.1;
	double searchRadius=scaleFactor*dx;
	CartesianComponent axis = UTILITIES::X;
	Array<int> bcIdsFixed = UTILITIES::getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,gridData.myX,gridData.numPoints,searchRadius,xMin);
	std::sort(bcIdsFixed.get(),bcIdsFixed.end());
	Array<int> bcIdsApplied = UTILITIES::getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,gridData.myX,gridData.numPoints,searchRadius,xMax);
	std::sort(bcIdsApplied.get(),bcIdsApplied.end());

	/**
	 * Create boundary conditions spec
	 */
	vector<shared_ptr<StageComponentDirichletBc> > bcs(2);
	ComponentDirichletBcSpec fixedSpec = ComponentDirichletBcSpec::getAllComponents(bcIdsFixed);
	StageFunction constStageFunction(0.0,0.0);
	shared_ptr<StageComponentDirichletBc> bcFixed(new StageComponentDirichletBc(fixedSpec,constStageFunction));
	bcs[0] = bcFixed;
	std::vector< DirichletBcSpec::ComponentLabel > c(1);
	c[0] = DirichletBcSpec::X;
	ComponentDirichletBcSpec appliedSpec(c,bcIdsApplied);
	StageFunction dispStageFunction(1.0e-3,1.0e-3);
	shared_ptr<StageComponentDirichletBc> bcApplied(new StageComponentDirichletBc(appliedSpec,dispStageFunction));
	bcs[1] = bcApplied;
	Field<char> bcMaskFieldOwned(BC_MASK,gridData.numPoints);
	bcMaskFieldOwned.set(0);
	for(int b=0;b<bcs.size();b++)
		bcs[b]->imprint_bc(bcMaskFieldOwned);
	Field<char> bcMaskFieldOverlap = PDNEIGH::createOverlapField(row_matrix_list,bcMaskFieldOwned);

	/*
	 * Create Jacobian -- note that SCOPE of jacobian is associated with the PimpOperator "op"
	 */
	Field<double> uOwnedField(DISPL3D,gridData.numPoints);
	uOwnedField.set(0.0);
	for(int b=0;b<bcs.size();b++)
		bcs[b]->applyKinematics(1.0,uOwnedField);
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op.getJacobian(uOwnedField);

	/*
	 * Create graph
	 */
	shared_ptr<Epetra_CrsGraph> graphPtr = getGraph(jacobian);

	/*
	 * Create Epetra_RowMatrix
	 */
//	shared_ptr<Epetra_RowMatrix> mPtr = getOperator(bcs,graphPtr,jacobian);
	shared_ptr<Epetra_RowMatrix> mPtr = getOperator_NEWER(bcMaskFieldOverlap,graphPtr,jacobian);

	/*
	 * TODO
	 * Investigate SIGN of applied loading
	 * 1) Compute internal force with displacement vector that has kinematics applied
	 * 2) Apply kinematics to this vector so that solution properly includes the applied kinematics
	 */
	Field_NS::FieldSpec fNSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fN");
	Field_NS::Field<double> fNOwnedField(fNSpec,gridData.numPoints);
	fNOwnedField.set(0.0);
	Field_NS::FieldSpec deltaSpec(FieldSpec::VELOCITY,FieldSpec::VECTOR3D,"delta");
	Field_NS::Field<double> delta(deltaSpec,gridData.numPoints);
	delta.set(0.0);
	op.computeInternalForce(uOwnedField,fNOwnedField);
	for(int b=0;b<bcs.size();b++){
		bcs[b]->applyKinematics(1.0,fNOwnedField);
		bcs[b]->applyKinematics(1.0,delta);
	}


	Epetra_LinearProblem linProblem;
	linProblem.SetOperator(mPtr.get());
	linProblem.AssertSymmetric();

	const Epetra_BlockMap& rangeMap  = mPtr->OperatorRangeMap();
	const Epetra_BlockMap& domainMap  = mPtr->OperatorDomainMap();

	double *f = fNOwnedField.get();
	Epetra_Vector rhs(View,rangeMap,f);
	Epetra_Vector lhs(View,domainMap,uOwnedField.get());

	linProblem.SetRHS(&rhs);
	linProblem.SetLHS(&lhs);
//	BOOST_CHECK(0==linProblem.CheckInput());

	AztecOO solver(linProblem);
	solver.SetAztecOption(AZ_precond, AZ_Jacobi);
//	BOOST_CHECK(0==solver.CheckInput());
	solver.Iterate(500,1e-6);
	/*
	 * Write problem set up parameters to file
	 */

	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(gridData.myX,gridData.numPoints);

	/*
	 * Add processor rank to output
	 */
	Field_NS::Field<int> myRank(Field_NS::PROC_NUM,gridData.numPoints);
	myRank.set(comm->MyPID());

	PdVTK::writeField(grid,fNOwnedField);
	PdVTK::writeField(grid,uOwnedField);
	PdVTK::writeField(grid,delta);
	PdVTK::writeField<int>(grid,myRank);
	PdVTK::writeField<char>(grid,bcMaskFieldOwned);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("utCrackOpeningDemo_npX.pvtu", comm->NumProc(), comm->MyPID());
	PdVTK::write(writer,grid);

}


int main
(
		int argc,
		char* argv[]
)
{
	// Initialize MPI and timer
	PdutMpiFixture myMpi = PdutMpiFixture(argc,argv);

	// These are static (file scope) variables
	myRank = myMpi.rank;
	numProcs = myMpi.numProcs;

	crackOpeningDemo();

}
