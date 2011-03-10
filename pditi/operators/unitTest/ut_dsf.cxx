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
#include "VectorUtils.h"

using namespace boost::unit_test;
using std::size_t;
using std::tr1::shared_ptr;
using namespace PdVTK;
using namespace Field_NS;
using namespace PdMaterialUtilities;
using namespace PdQuickGrid;
using std::pair;


template<class T> struct ArrayDeleter{
	void operator()(T* d) {
		delete [] d;
	}
};

template<class T>
class Array {

public:
	struct Deleter;
	friend struct Deleter;
	struct Deleter{
		void operator()(T* d) {
			delete [] d;
		}
	};

public:
	Array(size_t length) : aPtr(new T[length], Deleter()), size(length) {}
	size_t get_size() const { return size; }
	T* get() { return aPtr.get(); }
	const T* get() const { return aPtr.get(); }
	void set(T value) {
		T* s = aPtr.get();
		const T* e = s + size;
		for(; s!=e; s++)
			*s = value;
	}

private:
	shared_ptr<T> aPtr;
	size_t size;

};

void probe_shear
(
		PURE_SHEAR mode,
		Array<int> neighborhoodPtr,
		Array<double> X,
		Array<double> xPtr,
		Array<double> Y,
		Array<double> yPtr,
		Array<double> volPtr,
		double horizon,
		double gamma,
		double m_code
);

Array<double> create_spherical_neighborhood(const double center[], double radius, size_t numX, size_t numY, size_t numZ) {
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
	const int nx = numX;
	const int ny = numY;
	const int nz = numZ;
	size_t N = numX * numY * numZ;
	const double xStart  = -radius;
	const double xLength = 2 * radius;
	const double yStart  = -radius;
	const double yLength = 2 * radius;
	const double zStart  = -radius;
	const double zLength = 2 * radius;
	const PdQPointSet1d xSpec(nx,xStart,xLength);
	const PdQPointSet1d ySpec(ny,yStart,yLength);
	const PdQPointSet1d zSpec(nz,zStart,zLength);
	shared_ptr<double> xPtr=getDiscretization(xSpec,ySpec,zSpec);
	shared_ptr<bool> flagPtr(new bool[N],ArrayDeleter<bool>());
	size_t num_points = 0;
	{
		{
			double *X = xPtr.get();
			bool *flags = flagPtr.get();
			for(int p=0;p<N;p++, X+=3, flags++){
				*flags=false;
				VectorUtilsNS::Vector3D u(X);
				if(u.norm()<=radius){
					*flags=true;
					num_points+=1;
					*(X+0) += center[0];
					*(X+1) += center[1];
					*(X+2) += center[2];
				}
			}
		}
	}

	/*
	 * Now populate final array
	 */
	Array<double> a(3*num_points);
	{
		double *X = xPtr.get();
		double *X_Final = a.get();
		bool *flags = flagPtr.get();
		for(int p=0; p<N; p++, X+=3, flags++){
			if(*flags){
				*(X_Final+0) = *(X+0);
				*(X_Final+1) = *(X+1);
				*(X_Final+2) = *(X+2);
				X_Final+=3;
			}
		}
	}

	/*
	 * Write file for debugging and visualizing sphere of points
	 */
	int numProcs=1;
	int myRank=0;
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(a.get(), num_points);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer= PdVTK::getWriter("ut_dsf.pvtu", numProcs, myRank, PdVTK::vtkBINARY);
	PdVTK::write(writer,grid);

	return a;
}


void dsf_probe() {
	/*
	 * X is the center of the sphere
	 */
	Array<double> X(3); X.set(0.0);
	/*
	 * Y = X since we are fixing the center of the sphere
	 */
	Array<double> Y(3); Y.set(0.0);

	/*
	 * Radius of sphere
	 */
	const double radius = 1.0;
	const size_t num_points_along_axis(100);
	size_t numX(num_points_along_axis), numY(num_points_along_axis), numZ(num_points_along_axis);
	Array<double> xPtr = create_spherical_neighborhood(X.get(),radius,numX,numY,numZ);
	size_t num_points = xPtr.get_size()/3;
	Array<double> volPtr(num_points);
	double volSphere = 4.0 * M_PI * radius * radius * radius / 3.0;
	double avgPointVolume = volSphere / num_points;
	volPtr.set(avgPointVolume);

	/*
	 * In this case, the neighborhood is simply the list of point ids -- identity map
	 * NOTE: but need to put number of neighbors at the start and length of list
	 * is 'num_points+1'
	 */
	Array<int> neighborhoodPtr(num_points+1);
	{
		int *i = neighborhoodPtr.get();
		*i = num_points; i++;
		int *end = i + num_points;
		for(int j=0;i!=end;j++,i++)
			*i = j;
	}


	double m_analytical = 4.0 * M_PI * pow(radius,5) / 5.0;
	double m_code = computeWeightedVolume(X.get(),xPtr.get(),volPtr.get(),neighborhoodPtr.get());
	double rel_diff = abs(m_analytical-m_code)/m_analytical;
	double tolerance=1.0e-3;
	BOOST_CHECK_SMALL(rel_diff,tolerance);
//	std::cout << "ut_dsf::analytical value for weighted volume on sphere = " << m_analytical << std::endl;
//	std::cout << "ut_dsf::code computed weighted volume on sphere = " << m_code << std::endl;

	double gamma = 1.0e-6;
	double horizon = radius;
	Array<double> yPtr(3*num_points);

	/*
	 * PROBE XY
	 */
	probe_shear(XY,neighborhoodPtr,X,xPtr,Y,yPtr,volPtr,horizon,gamma,m_code);
	/*
	 * PROBE XZ
	 */
	probe_shear(XZ,neighborhoodPtr,X,xPtr,Y,yPtr,volPtr,horizon,gamma,m_code);

	/*
	 * PROBE YZ
	 */
	probe_shear(YZ,neighborhoodPtr,X,xPtr,Y,yPtr,volPtr,horizon,gamma,m_code);



}

void probe_shear
(
	PURE_SHEAR mode,
	Array<int> neighborhoodPtr,
	Array<double> X,
	Array<double> xPtr,
	Array<double> Y,
	Array<double> yPtr,
	Array<double> volPtr,
	double horizon,
	double gamma,
	double m_code
)
{
	/*
	 * This is the reference value for ed_squared
	 */
	double reference = 4.0 * M_PI * gamma * gamma * pow(horizon,5) / 75.0;
	/*
	 * NOTE: X is center of sphere and there no displacement at this point
	 * therefore, Y=X
	 */
	set_pure_shear(neighborhoodPtr.get(),X.get(),xPtr.get(),yPtr.get(),mode,gamma);
	double theta = computeDilatation(neighborhoodPtr.get(),X.get(),xPtr.get(),X.get(),yPtr.get(),volPtr.get(),m_code);
	double tolerance=1.0e-12;
	BOOST_CHECK_SMALL(theta,tolerance);

	/*
	 * compute shear correction factor
	 */
	double ed_squared = compute_norm_2_deviatoric_extension(neighborhoodPtr.get(),X.get(),xPtr.get(),Y.get(),yPtr.get(),volPtr.get());
	double dsf = reference/ed_squared;
	//	std::cout << "ut_dsf::probe_shear MODE = " << mode << std::endl;
	//	std::cout << "ut_dsf::probe_shear computed dilatation in pure shear = " << theta << std::endl;
	/*
	 * For this nearly perfect 'sphere', the shear correction factor should be very close to '1.0'
	 */
	double rel_diff = abs(1.0-dsf);
	tolerance=1.0e-3;
	BOOST_CHECK_SMALL(rel_diff,tolerance);

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
