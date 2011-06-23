/*
 * ut_OverlapDistributor.cxx
 *
 *  Created on: Mar 30, 2011
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "Sortable.h"
#include "Array.h"
#include "quick_grid/QuickGrid.h"
#include "../NeighborhoodList.h"
#include "../OverlapDistributor.h"
#include "../BondFilter.h"
#include "../PdZoltan.h"
#include "vtk/Field.h"
#include "../../operators/DirichletBcSpec.h"
#include "../../operators/StageFunction.h"
#include "../../operators/ComponentDirichletBc.h"
#include "../../operators/ComponentDirichletBcSpec.h"

#include "PdutMpiFixture.h"
#include <iostream>
#include <sstream>
#include <vector>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


using UTILITIES::CartesianComponent;
using UTILITIES::Array;
using namespace PdBondFilter;
using PdITI::STAGE::ComponentDirichletBc;
using PdITI::STAGE::DirichletBc;
using PdITI::STAGE::DirichletBcSpec;
using PdITI::STAGE::ComponentDirichletBcSpec;
using PdITI::STAGE::StageFunction;
using namespace PdITI::Field_NS;
using namespace Pdut;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::vector;
using std::cout;
using std::endl;
using std::stringstream;

const int vectorNDF=3;
static size_t myRank;
static size_t numProcs;

/*
 * This should be even so that the crack plane lies between to rows of points
 */
const size_t nx = 4;
const size_t ny = 4;
const size_t nz = 1;
const double xStart = -2.5;
const double xLength = 5.0;
const double yStart = -2.5;
const double yLength = 5.0;
const double zStart = -0.5;
const double zLength = 1.0;
static double xMax=xStart+xLength;
static double xMin=xStart;
const size_t numCells = nx*ny*nz;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const double dx = xSpec.getCellSize();
const double dy = ySpec.getCellSize();
const double dz = zSpec.getCellSize();
const double _cellVolume = dx*dy*dz;


/*
 * Horizon
 */
const double horizon=1.1*sqrt( (3.0*dx)*(3.0*dx) );

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

void dirichletBc() {
    /*
     * Get mesh and decomposition
     */
    QUICKGRID::QuickGridData gridData = getGrid();
    /*
     * Communicator
     */
    shared_ptr<Epetra_Comm> comm = shared_ptr<Epetra_Comm>(new Epetra_MpiComm(MPI_COMM_WORLD));

    shared_ptr<BondFilter> filterPtr = shared_ptr<BondFilter>(new BondFilterDefault());
    PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);
    PDNEIGH::NeighborhoodList row_matrix_list = list.cloneAndShare(2.0 * list.get_horizon());
    const Epetra_BlockMap& ownedMap = *row_matrix_list.getOwnedMap(1);
    const Epetra_BlockMap& overlapMap = *row_matrix_list.getOverlapMap(1);

    /*
     * Get points for bc's
     *
     * Note that we are looking for a discrete number of points at end;
     * Set the scale factor to just larger than an integer where
     * the integer corresponds with the number of points to be included in
     * the boundary conditions
     */
    double scaleFactor=1.1;
    double searchRadius=scaleFactor*dx;
    CartesianComponent axis = UTILITIES::X;
    Array<int> bcIdsFixed = UTILITIES::getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,gridData.myX,gridData.numPoints,searchRadius,xMin);
    Array<int> bcIdsApplied = UTILITIES::getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,gridData.myX,gridData.numPoints,searchRadius,xMax);

    /**
     * Create boundary conditions spec
     */
    vector<shared_ptr<ComponentDirichletBc> > bcs(2);
    ComponentDirichletBcSpec fixedSpec = ComponentDirichletBcSpec::getAllComponents(bcIdsFixed);
    StageFunction constStageFunction(0.0,0.0);
    shared_ptr<ComponentDirichletBc> bcFixed(new ComponentDirichletBc(fixedSpec,constStageFunction));
    bcs[0] = bcFixed;
    std::vector< DirichletBcSpec::ComponentLabel > c(1);
    c[0] = DirichletBcSpec::X;
    ComponentDirichletBcSpec appliedSpec(c,bcIdsApplied);
    StageFunction dispStageFunction(1.0e-3,1.0e-3);
    shared_ptr<ComponentDirichletBc> bcApplied(new ComponentDirichletBc(appliedSpec,dispStageFunction));
    bcs[1] = bcApplied;
    Field<char> bcMaskFieldOwned(PdITI::BC_MASK,gridData.numPoints);
    bcMaskFieldOwned.set(0);
    for(int b=0;b<bcs.size();b++)
        bcs[b]->imprint_bc(bcMaskFieldOwned);
    Field<char> bcMaskFieldOverlap = PDNEIGH::createOverlapField(row_matrix_list,bcMaskFieldOwned);

    BOOST_CHECK(bcMaskFieldOverlap.get_num_points()==16);

    /*
     * All GIDs are in the overlap map on every processor
     */
    int GIDS[] =  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    char bcMask[]={7,0,0,1,7,0,0,1,7,0, 0, 1, 7, 0, 0, 1};
    stringstream streamGIDs;
    for(int i=0;i<bcMaskFieldOverlap.get_num_points(); i++){
        int lid = overlapMap.LID(GIDS[i]);
        int mask = bcMaskFieldOverlap[lid];
        streamGIDs <<  GIDS[i] << ", " << lid << ", " << mask  << "\n";
        BOOST_CHECK(bcMask[i]==mask);
    }
    cout << streamGIDs.str() << endl;
    if(0==myRank){

        /*
         * Ids at fixed end
         */
        BOOST_CHECK(bcIdsFixed.get_size()==0);

        /*
         * Ids at applied end
         */
        BOOST_CHECK(bcIdsApplied.get_size()==4);
        int bcIdsApplied_answer[]={3,7,11,15};
        for(int i=0;i<4;i++){
            BOOST_CHECK(bcIdsApplied_answer[i]==ownedMap.GID(bcIdsApplied[i]));
        }
    }
    if(1==myRank){
        /*
         * Ids at fixed end
         */
        BOOST_CHECK(bcIdsFixed.get_size()==2);
        int bcIdsFixed_answer[]={0,4};
        for(int i=0;i<2;i++)
            BOOST_CHECK(bcIdsFixed_answer[i]==ownedMap.GID(bcIdsFixed[i]));

        /*
         * Ids at applied end
         */
        BOOST_CHECK(bcIdsApplied.get_size()==0);
    }
    if(2==myRank){
        /*
         * Ids at fixed end
         */
        BOOST_CHECK(bcIdsFixed.get_size()==2);
        int bcIdsFixed_answer[]={12,8};
        for(int i=0;i<2;i++){
            BOOST_CHECK(bcIdsFixed_answer[i]==ownedMap.GID(bcIdsFixed[i]));
        }

        /*
         * Ids at applied end
         */
        BOOST_CHECK(bcIdsApplied.get_size()==0);

    }



}

bool init_unit_test_suite()
{
    // Add a suite for each processor in the test
    bool success=true;

    test_suite* proc = BOOST_TEST_SUITE( "dirichletBc" );
    proc->add(BOOST_TEST_CASE( &dirichletBc ));
    framework::master_test_suite().add( proc );
    return success;

}

bool init_unit_test()
{
    init_unit_test_suite();
    return true;
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
    /**
     * This test only make sense for numProcs == 3
     */
    if(3 != numProcs){
        std::cerr << "Unit test runtime ERROR: ut_OverlapDistributor only makes sense on 4 processors" << std::endl;
        std::cerr << "\t Re-run unit test $mpiexec -np 3 ./ut_OverlapDistributor." << std::endl;
        myMpi.PdutMpiFixture::~PdutMpiFixture();
        std::exit(-1);
    }

    // Initialize UTF
    return unit_test_main( init_unit_test, argc, argv );
}
