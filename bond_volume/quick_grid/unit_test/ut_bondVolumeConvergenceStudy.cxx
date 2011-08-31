/*
 * ut_bondVolumeConvergenceStudy.cxx
 *
 *  Created on: Aug 2, 2011
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <cmath>
#include <tr1/memory>
#include <iostream>
#include <fstream>
#include <sstream>


#include "mesh_output/vtk/PdVTK.h"
#include "mesh_output/vtk/Field.h"
#include "mesh_input/quick_grid/QuickGrid.h"
#include "bond_volume/quick_grid/calculators.h"
#include "ordinary_utilities.h"
#include "pdneigh/NeighborhoodList.h"
#include "pdneigh/PdZoltan.h"
#include "pdneigh/BondFilter.h"
#include "utilities/Array.h"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

using namespace boost::unit_test;
using std::size_t;
using std::tr1::shared_ptr;
using namespace PdVTK;
using UTILITIES::Array;
using UTILITIES::Vector3D;
using std::pair;



static int nx;
static int ny;
static int nz;
const size_t numProcs=1;
const size_t myRank=0;


void probe_shear
(
		MATERIAL_EVALUATION::PURE_SHEAR mode,
		Array<int> neighborhoodPtr,
		Array<double> X,
		Array<double> xPtr,
		Array<double> Y,
		Array<double> yPtr,
		Array<double> bondVolume,
		double horizon,
		double gamma,
		double m_code
);

void dsf_probe(const std::string& json_filename);

void set_static_data(const std::string& json_filename)
{
	// Create an empty property tree object
	using boost::property_tree::ptree;
	ptree pt;

	try {
		read_json(json_filename, pt);
	} catch(std::exception& e){
		std::cerr << e.what();
		std::exit(1);
	}

	/*
	 * Get Discretization
	 */
	ptree discretization_tree=pt.find("Discretization")->second;
	std::string path=discretization_tree.get<std::string>("Type");

	if("QuickGrid.TensorProduct3DMeshGenerator"==path){

		nx = pt.get<int>(path+".Number Points X");
		ny = pt.get<int>(path+".Number Points Y");
		nz = pt.get<int>(path+".Number Points Z");

	} else {
		std::string s;
		s = "Error-->ut_bondVolumeConvergenceStudy\n";
		s += "\tTest only works for Discretization.Type==QuickGrid.TensorProduct3DMeshGenerator\n";
		throw std::runtime_error(s);
	}

}


void write_table_1_header(const std::string& output_tex_table){
	std::stringstream table_out;

	table_out << "\\begin{table}[ht]" << "\n";
	table_out << "\\centering" << "\n";
	table_out << "\\bigskip" << "\n";
	table_out << "\\begin{tabular}{|c|c|c|c|}" << "\n";
	table_out << "\\hline" << "\n";
	table_out << "$n$ "
			    << "& $\\frac{|m-m_n|}{m}$ "
			    << "& $\\frac{\\Vert e^d\\Vert^2-\\Vert e^d_n\\Vert^2}{\\Vert e^d\\Vert^2}$ "
			    << "& $\\frac{\\Vert e^d\\Vert^2}{\\Vert e^d_n\\Vert^2}$ \\\\" << "\n";
	table_out << "\\hline" << "\n";


	std::ofstream file_stream;
	file_stream.open(output_tex_table.c_str(),ios::app|ios::out);

	file_stream << table_out.str();
	file_stream.close();

}

void close_table_1(const std::string& output_tex_table) {
	std::stringstream table_out;

	table_out << "\\hline" << "\n";
	table_out << "\\end{tabular}" << "\n";
	table_out << "\\end{table}" << "\n";
	std::ofstream file_stream;
	file_stream.open(output_tex_table.c_str(),ios::app|ios::out);
	file_stream << table_out.str();
	file_stream.close();
}


QUICKGRID::QuickGridData getGrid(const string& _json_filename) {
	shared_ptr<QUICKGRID::QuickGridMeshGenerationIterator> g;
	g = QUICKGRID::getMeshGenerator(numProcs,_json_filename);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, *g);

	// This load-balances
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
	return decomp;
}

void ut_bondVolumeConvergenceStudy_n3() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=3.json";
	set_static_data(file);
	dsf_probe(file);
}

void ut_bondVolumeConvergenceStudy_n5() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=5.json";
	set_static_data(file);
	dsf_probe(file);
}

void ut_bondVolumeConvergenceStudy_n7() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=7.json";
	set_static_data(file);
	dsf_probe(file);
}

void ut_bondVolumeConvergenceStudy_n9() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=9.json";
	set_static_data(file);
	dsf_probe(file);
}

void ut_bondVolumeConvergenceStudy_n11() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=11.json";
	set_static_data(file);
	dsf_probe(file);
}

void ut_bondVolumeConvergenceStudy_n13() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=13.json";
	set_static_data(file);
	dsf_probe(file);
}

void ut_bondVolumeConvergenceStudy_n17() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=17.json";
	set_static_data(file);
	dsf_probe(file);
}

void ut_bondVolumeConvergenceStudy_n33() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=33.json";
	set_static_data(file);
	dsf_probe(file);
}



/*
 * Dave's Influence Function
 * "x < 0.5 ? 1.0 : -4.0*x*x + 4.0*x"
 */

void dsf_probe(const std::string& json_filename) {
//	const string json_filename="./input_files/ut_dsf.json";
	shared_ptr<BOND_VOLUME::QUICKGRID::Bond_Volume_Calculator> c;
	c = BOND_VOLUME::QUICKGRID::get_Bond_Volume_Calculator(json_filename);
	double horizon = c->get_horizon();

	QUICKGRID::QuickGridData gridData = getGrid(json_filename);

	// This load-balances
	gridData = PDNEIGH::getLoadBalancedDiscretization(gridData);

	/*
	 * Create neighborhood with an enlarged horizon
	 */
	double cell_diagonal=c->get_cell_diagonal();
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon+cell_diagonal);


	/*
	 * Unit test looks exclusively at the cell at the center of cube;
	 * This cell ID depends upon nx, ny, nz
	 *
	 * MESH INPUT MUST HAVE EVEN NUMBER OF CELLS
	 */
	BOOST_CHECK(0==(nx+1)%2);
	/*
	 * mesh must have equal number of cells along each axis
	 */
	BOOST_CHECK(nx==ny);
	BOOST_CHECK(nx==nz);

	// coordinate indices of center cell
	size_t ic = (nx -1)/2;
	size_t jc = ic;
	size_t kc = ic;
	size_t gId = nx * ny * kc + nx * jc + ic;
//	std::cout << "ut_dsf::center cell gID = " << gId << std::endl;

	/**
	 * WARNING: note following ASSUMPTION -- gId == local id
	 * CAUTION: this test only makes sense in 'serial' -- local id
	 * and gId are not the same in parallel.
	 */
	int num_neigh = list.get_num_neigh(gId);
//	std::cout << "ut_dsf::center cell num_neigh = " << num_neigh << std::endl;
	Array<int> neighborhoodPtr(1+num_neigh);
	{
		/*
		 * copy neighborhood list for center point over to array
		 */
		const int *neighborhood = list.get_neighborhood(gId);
		BOOST_CHECK(num_neigh == *neighborhood);
		for(size_t j=0;j<num_neigh+1;j++,neighborhood++)
			neighborhoodPtr[j]=*neighborhood;
	}
	Array<double> xPtr(list.get_num_owned_points()*3,list.get_owned_x());
//	shared_ptr<double> xPtr = list.get_owned_x();


	/*
	 * expectation is that cell center is at origin
	 */
//	std::cout << "ut_dsf::cell center coordinate X = " << *(xPtr.get()+3*gId) << std::endl;
//	std::cout << "ut_dsf::cell center coordinate Y = " << *(xPtr.get()+3*gId + 1) << std::endl;
//	std::cout << "ut_dsf::cell center coordinate Z = " << *(xPtr.get()+3*gId + 2) << std::endl;
	BOOST_CHECK_SMALL(xPtr[3*gId+0],1.0e-15);
	BOOST_CHECK_SMALL(xPtr[3*gId+1],1.0e-15);
	BOOST_CHECK_SMALL(xPtr[3*gId+2],1.0e-15);
	/*
	 * X is the center of the sphere
	 */
	Array<double> X(3); X.set(0.0);
	/*
	 * Y = X since we are fixing the center of the sphere
	 */
	Array<double> Y(3); Y.set(0.0);


	/*
	 * only computing bondVolume on center point 'X'
	 */
	Array<double> bondVolume(num_neigh);
	bondVolume.set(0.0);
	BOND_VOLUME::QUICKGRID::compute_bond_volume(X.get(),neighborhoodPtr.get(),xPtr.get(),bondVolume.get(),c.get());

	double m_analytical = 4.0 * M_PI * pow(horizon,5) / 5.0;
	double m_code = MATERIAL_EVALUATION::WITH_BOND_VOLUME::computeWeightedVolume(X.get(),xPtr.get(),bondVolume.get(),neighborhoodPtr.get());
	double rel_diff = std::abs(m_analytical-m_code)/m_analytical;
	std::cout << std::scientific;
	std::cout.precision(3);
	std::cout << "ut_dsf::analytical value for weighted volume on sphere = " << m_analytical << std::endl;
	std::cout << "ut_dsf::code computed weighted volume on sphere = " << m_code << std::endl;
	std::cout << "ut_dsf::% relative error weighted volume = " << 100*rel_diff << std::endl;

	double gamma = 1.0e-6;
	Array<double> yPtr(3*list.get_num_owned_points());

	/*
	 * PROBE XY
	 */
	probe_shear(MATERIAL_EVALUATION::XY,neighborhoodPtr,X,xPtr,Y,yPtr,bondVolume,horizon,gamma,m_code);
	/*
	 * PROBE YZ
	 */
//	probe_shear(MATERIAL_EVALUATION::YZ,neighborhoodPtr,X,xPtr,Y,yPtr,bondVolume,horizon,gamma,m_code);
	/*
	 * PROBE ZX
	 */
//	probe_shear(MATERIAL_EVALUATION::ZX,neighborhoodPtr,X,xPtr,Y,yPtr,bondVolume,horizon,gamma,m_code);


}

void probe_shear
(
	MATERIAL_EVALUATION::PURE_SHEAR mode,
	Array<int> neighborhoodPtr,
	Array<double> X,
	Array<double> xPtr,
	Array<double> Y,
	Array<double> yPtr,
	Array<double> bondVolume,
	double horizon,
	double gamma,
	double m_code
)
{

	/*
	 * This is the reference value for the weighted volume
	 */
	double m_analytical = 4.0 * M_PI * pow(horizon,5) / 5.0;
	double m_err = std::fabs(m_analytical-m_code)/m_analytical;
	/*
	 * NOTE: X is center of sphere and there no displacement at this point
	 * therefore, Y=X
	 */
	MATERIAL_EVALUATION::set_pure_shear(neighborhoodPtr.get(),X.get(),xPtr.get(),yPtr.get(),mode,gamma);
	double theta = MATERIAL_EVALUATION::WITH_BOND_VOLUME::computeDilatation(neighborhoodPtr.get(),X.get(),xPtr.get(),X.get(),yPtr.get(),bondVolume.get(),m_code);
//	std::cout << "ut_bondVolumeConvergenceStudy::probe_shear dilatation = " << theta << std::endl;
	double tolerance=1.0e-12;
	BOOST_CHECK_SMALL(theta,tolerance);

	/*
	 * compute shear correction factor
	 */
	/*
	 * This is the reference value for ed_squared
	 */
	double reference = 4.0 * M_PI * gamma * gamma * pow(horizon,5) / 75.0;
	double ed2 = MATERIAL_EVALUATION::WITH_BOND_VOLUME::compute_norm_2_deviatoric_extension(neighborhoodPtr.get(),X.get(),xPtr.get(),Y.get(),yPtr.get(),bondVolume.get(),m_code);
	double dsf = reference/ed2;
	double ed_err = fabs(reference-ed2)/reference;
	std::cout << "ut_dsf::probe_shear MODE = " << mode << std::endl;
	std::cout << "ut_dsf::ed = " << reference << std::endl;
	cout.precision(2);
	std::cout << std::scientific << "ut_dsf::probe_shear computed % ed_err in pure shear = " << 100*ed_err << std::endl;
	std::cout << "ut_dsf::probe_shear computed dsf in pure shear = " << dsf << std::endl;
	/*
	 * For this nearly perfect 'sphere', the shear correction factor should be very close to '1.0'
	 */

	std::stringstream table_1_out;
	table_1_out << nx << " & ";
	table_1_out.precision(4);
	table_1_out << m_err*100 << "\\% & ";
	table_1_out << ed_err*100 << "\\% & ";
	table_1_out.precision(3);
	table_1_out << dsf << " \\\\ \n";

	std::ofstream file_stream;
	file_stream.open("table_1.tex",ios::app|ios::out);
	file_stream << table_1_out.str();
	file_stream.close();

	/*
	 * write raw data
	 */
	file_stream.open("ut_bondVolumeConvergenceStudy.dat",ios::app|ios::out);
	file_stream << nx << " ";
	file_stream << std::scientific;
	file_stream.precision(15);
	file_stream << 2.0*horizon/nx << " ";
	file_stream << m_code << " ";
	file_stream << ed2 << "\n";
	file_stream.close();


}


bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "ut_bondVolumeConvergenceStudy" );
	proc->add(BOOST_TEST_CASE( &ut_bondVolumeConvergenceStudy_n3 ));
	proc->add(BOOST_TEST_CASE( &ut_bondVolumeConvergenceStudy_n5 ));
	proc->add(BOOST_TEST_CASE( &ut_bondVolumeConvergenceStudy_n7 ));
	proc->add(BOOST_TEST_CASE( &ut_bondVolumeConvergenceStudy_n9 ));
	proc->add(BOOST_TEST_CASE( &ut_bondVolumeConvergenceStudy_n11 ));
	proc->add(BOOST_TEST_CASE( &ut_bondVolumeConvergenceStudy_n13 ));
	proc->add(BOOST_TEST_CASE( &ut_bondVolumeConvergenceStudy_n17 ));
//	proc->add(BOOST_TEST_CASE( &ut_bondVolumeConvergenceStudy_n33 ));
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


	write_table_1_header("table_1.tex");
//	write_table_2_header("table_2.tex");

	// Initialize UTF
	int flag = unit_test_main( init_unit_test, argc, argv );

	close_table_1("table_1.tex");
	return flag;
}
