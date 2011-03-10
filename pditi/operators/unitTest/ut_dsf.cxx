/*
 * ut_dsf.cxx
 *
 *  Created on: Mar 9, 2011
 *      Author: jamitch
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <math.h>
#include <tr1/memory>
#include "PdVTK.h"
#include "Field.h"
#include "PdMaterialUtilities.h"
#include "../PdITI_Utilities.h"
#include "PdQuickGrid.h"

using namespace boost::unit_test;
using std::size_t;
using std::tr1::shared_ptr;
using namespace PdVTK;
using namespace Field_NS;
using namespace PdMaterialUtilities;
using namespace PdQuickGrid;

template<class T> struct ArrayDeleter{
	void operator()(T* d) {
		delete [] d;
	}
};

shared_ptr<double> create_spherical_neighborhood(double center[], double radius, size_t numR, size_t numTheta, size_t numPhi) {
	/*
	 * Random coordinates within a sphere of 'radius' and centered at 'center'
	 */
	size_t N = numR * numTheta * numPhi;
	shared_ptr<double> xPtr(new double[3*N],ArrayDeleter<double>());
	{
		/*
		 * Initialize random number generator
		 */
		srand ( time(NULL) );


		double pi = M_PI;
		/*
		 * Spherical coordinates
		 * x = r * sin(theta) * cos(phi)
		 * y = r * sin(theta) * sin(phi)
		 * z = r * cos(theta)
		 *
		 * theta in [0,Pi]
		 * phi   in [0,2 Pi)
		 */
//		PdQPointSet1d specR(numR,0,radius);
//		PdQPointSet1d specTheta(numTheta,0,M_PI);
//		PdQPointSet1d specPhi(numPhi,0,2.0*M_PI);
//		shared_ptr<double> thetaPtr = PdQuickGrid::getDiscretization(specTheta);
//		shared_ptr<double> phiPtr = PdQuickGrid::getDiscretization(specPhi);
//		shared_ptr<double> rPtr = PdQuickGrid::getDiscretization(specR);

		double *X = xPtr.get();
		for(int nR=0;nR<numR;nR++){
//			double r     = *(rPtr.get()+nR);
			double r     = radius * (rand()%numR)/numR;
			for(int nTheta=0;nTheta<numTheta;nTheta++){
//				double theta = *(thetaPtr.get()+nTheta);
				double theta     = M_PI * (rand()%numTheta)/numTheta;
				for(int nPhi=0;nPhi<numPhi;nPhi++){
//					double phi   = *(phiPtr.get()+nPhi);
					double phi     = 2.0 * M_PI * (rand()%numPhi)/numPhi;
					*(X+0)= center[0] + r * sin(theta) * cos(phi);
					*(X+1)= center[1] + r * sin(theta) * sin(phi);
					*(X+2)= center[2] + r * cos(theta);
					X+=3;
				}

			}

		}
	}

	/*
	 * Write file for debugging and visualizing sphere of points
	 */
	int numProcs=1;
	int myRank=0;
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(xPtr.get(), N);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer= PdVTK::getWriter("ut_dsf.pvtu", numProcs, myRank, PdVTK::vtkBINARY);
	PdVTK::write(writer,grid);

	return xPtr;
}


void dsf_probe() {
	double x0[] = {0.0,0.0,0.0};
	double radius = 1.0;
	size_t numR(100), numTheta(100), numPhi(100);
	size_t num_points = numR * numTheta * numPhi;
	shared_ptr<double> xPtr = create_spherical_neighborhood(x0,radius,numR,numTheta,numPhi);
	shared_ptr<double> volPtr(new double[num_points],ArrayDeleter<double>());
	double volSphere = 4.0 * M_PI * radius * radius * radius / 3.0;
	double avgPointVolume = volSphere / num_points;
	PdITI::SET(volPtr.get(),volPtr.get()+num_points,avgPointVolume);
	std::cout << "ut_dsf::analytical value for volume of sphere = " << volSphere << std::endl;
	std::cout << "ut_dsf::volume point in sphere = " << avgPointVolume << std::endl;

	/*
	 * In this case, the neighborhood is simply the list of point ids -- identity map
	 * NOTE: but need to put number of neighbors at the start and length of list
	 * is 'num_points+1'
	 */
	shared_ptr<int> neighborhoodPtr(new int[num_points+1],ArrayDeleter<int>());
	{
		int *i = neighborhoodPtr.get();
		*i = num_points; i++;
		int *end = i + num_points;
		for(int j=0;i!=end;j++,i++)
			*i = j;

	}


	double m_analytical = 4.0 * M_PI * pow(radius,5) / 5.0;
	double m_code = computeWeightedVolume(x0,xPtr.get(),volPtr.get(),neighborhoodPtr.get());
	std::cout << "ut_dsf::analytical value for weighted volume on sphere = " << m_analytical << std::endl;
	std::cout << "ut_dsf::code computed weighted volume on sphere = " << m_code << std::endl;

	double gamma = 1.0e-6;
	shared_ptr<double> yPtr(new double[3*num_points],ArrayDeleter<double>());
//	set_pure_shear(neighborhoodPtr.get(),x0,xPtr.get(),yPtr.get(),XY,gamma);
//
//
//	double horizon = radius;
//	/*
//	 * there is no displacement at center point; therefore y0 = x0
//	 */
//	double y0[] = {x0[0],x0[1],x0[2]};
//	double dsf = probeShearModulusScaleFactor(neighborhoodPtr.get(),x0,xPtr.get(),y0,yPtr.get(),volPtr.get(),horizon,gamma);
//	std::cout << "ut_dsf::dsf = " << dsf << std::endl;
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "ut_dsf" );
	proc->add(BOOST_TEST_CASE( &dsf_probe ));
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
