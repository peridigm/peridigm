
#include <Shards_CellTopology.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_Cubature.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Intrepid_HGRAD_QUAD_C1_FEM.hpp>
#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>
#include <Teuchos_RCP.hpp>
#include <iostream>

using Teuchos::ArrayRCP;
using Teuchos::RCP;
using std::cout;
using std::endl;


typedef Intrepid::FieldContainer<double> FieldContainer;

typedef shards::CellTopology CellTopology;

typedef Intrepid::Cubature<double> Cubature;

typedef std::size_t size_t;
static size_t numFields = 4;
static size_t spaceDim = 2;

class MyData {

public:
	MyData(std::size_t numCells, std::size_t cubDegree, const CellTopology& cellType)
	: numCells(numCells), cubDegree(cubDegree), num_cubature_points(-1), cellType(cellType)
	{
		Intrepid::DefaultCubatureFactory<double> cubFactory;
		myCub = cubFactory.create(cellType, cubDegree);
		num_cubature_points = myCub->getNumPoints();

		/*
		 * basis at quadrature points
		 */
		basis_at_cub_points = RCP<FieldContainer>(new FieldContainer(numFields, num_cubature_points));
		grad_at_cub_points  = RCP<FieldContainer>(new FieldContainer(numFields, num_cubature_points, spaceDim));

		// Allocate and format arrays
		cub_points =   RCP<FieldContainer>(new FieldContainer           (num_cubature_points, spaceDim));
		cub_weights =  RCP<FieldContainer>(new FieldContainer           (num_cubature_points));

		// Compute basis functions and gradients at quadrature points
		myCub->getCubature(*cub_points, *cub_weights);
		quadBasis.getValues(*basis_at_cub_points, *cub_points, Intrepid::OPERATOR_VALUE);
		quadBasis.getValues(*grad_at_cub_points,  *cub_points, Intrepid::OPERATOR_GRAD);

	}

	~MyData() { cout << "~MyData()" << std::endl; }

private:
    std::size_t numCells, cubDegree, num_cubature_points;

    const CellTopology cellType;

    RCP<Cubature> myCub;

	// HGRAD basis on QUAD
	Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, FieldContainer > quadBasis;

    /*
     * F,P
     * numFields, num_cubature_points
     */
    RCP<FieldContainer> basis_at_cub_points;

    /*
     * F, P, D
     */
    RCP<FieldContainer> grad_at_cub_points;

    /*
     * P,D
     */
    RCP<FieldContainer> cub_points;

    /*
     * P
     */
    RCP<FieldContainer> cub_weights;


};

int main
(
		int argc,
		char* argv[]
)
{


	std::size_t numCells = 1;
	std::size_t cubDegree=2;
	shards::CellTopology cellType = shards::getCellTopologyData< shards::Quadrilateral<4> >();
	MyData myData(numCells,cubDegree,cellType);
}
