/*
 * utFIFE_2D.cxx
 *
 *  Created on: Oct 15, 2010
 *      Author: awesome
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <Shards_CellTopology.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_Cubature.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Intrepid_HGRAD_QUAD_C1_FEM.hpp>
#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>
#include "../FIFE_2D_IntrepidQuadrature.h"

using namespace boost::unit_test;
using namespace FIFE_2D_IntrepidQuadrature;

static size_t numFields = 4;
static size_t spaceDim = 2;

void FIFE_2D()
{
	shards::CellTopology cellType = shards::getCellTopologyData< shards::Quadrilateral<4> >();
	std::size_t numCells = 10000;
	std::size_t npe      = cellType.getNodeCount();
	std::size_t cubDegree=2;
	IntrepidData iData(numCells,cubDegree,cellType);

	// Allocate and format arrays
	FieldContainer cell_nodes(numCells, npe, spaceDim);

	for(size_t e=0;e<numCells;e++) {
		/*
		 * Create 'unit' cell
		 */
		cell_nodes(e,0,0)=-1.0; cell_nodes(e,0,1)=-1.0;
		cell_nodes(e,1,0)= 1.0; cell_nodes(e,1,1)=-1.0;
		cell_nodes(e,2,0)= 1.0; cell_nodes(e,2,1)= 1.0;
		cell_nodes(e,3,0)=-1.0; cell_nodes(e,3,1)= 1.0;

	}

	/*
	 * Mean magnetic field
	 */
	double B0[2] = {0.0,0.0};

	/*
	 * A3Np1 -- vector potential at step n+1 (C,F)
	 */
	FieldContainer A3Np1(numCells, numFields);
	double A_BOTTOM=0.0;
	double A_TOP   =1.0;
	for(size_t e=0;e<numCells;e++) {
		A3Np1(e,0)=A_BOTTOM; A3Np1(e,1)=A_BOTTOM;
		A3Np1(e,2)=A_TOP;    A3Np1(e,3)=A_TOP;
	}

	/*
	 * initialize intrepid data and compute ||B||
	 */
	iData.computeNormB(cell_nodes,A3Np1,B0);

    /*
     * Assign material properties
     */
    FieldContainer& sigmaGP = iData.getConductivity();
    FieldContainer& nuGP = iData.getReluctivity();
    FieldContainer& B = iData.getNormB();
    for(size_t e=0;e<numCells;e++)
    	for(size_t p=0;p<iData.getNumCubaturePoints();p++){
    		double b = B(e,p);
    		sigmaGP(e,p)=1.0;
    		nuGP(e,p)=1.0;
    	}


	/*
	 * Mass and stiffness matrices
	 */
	FieldContainer stiffness_matrices(numCells, numFields, numFields);
	FieldContainer mass_matrices(numCells, numFields, numFields);

	/*
	 * complete quadrature
	 */
	iData.finish_quadrature(stiffness_matrices,mass_matrices);

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utLinkIntrepid" );
	proc->add(BOOST_TEST_CASE( &FIFE_2D ));
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


	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
