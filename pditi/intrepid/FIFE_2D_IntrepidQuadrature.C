
//#ifdef SNL_MHD
//
//#ifdef TWO_D
//
//#ifdef HAVE_INTREPID

#include "FIFE_2D_IntrepidQuadrature.h"
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Intrepid_HGRAD_QUAD_C1_FEM.hpp>
#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

namespace FIFE_2D_IntrepidQuadrature {

using std::size_t;
static size_t numFields = 4;
static size_t spaceDim = 2;



RCP<FieldContainer> computeCellArea(const FieldContainer& n, const CellTopology& cellType) {

	size_t numCells = n.dimension(0);
	RCP<FieldContainer> area = RCP<FieldContainer>(new FieldContainer (numCells));
	FieldContainer a = *area;

	for (size_t e=0;e<numCells;e++){
		double x0=n(e,0,0);
		double x1=n(e,1,0);
		double x2=n(e,2,0);
		double x3=n(e,3,0);
		double y0=n(e,0,1);
		double y1=n(e,1,1);
		double y2=n(e,2,1);
		double y3=n(e,3,1);
		a(e) = 0.5 * ((x3-x1)*(y0-y2)+(x0-x2)*(y1-3));
	}

	return area;
}


RCP<FieldContainer> meanQuadratureCurlPhi(const FieldContainer& n, const CellTopology& cellType) {

	size_t numCells = n.dimension(0);
	RCP<FieldContainer> curlPhi = RCP<FieldContainer>(new FieldContainer(numCells,numFields,spaceDim));
	FieldContainer curl = *curlPhi;

	for (size_t e=0;e<numCells;e++){
		double x0=n(e,0,0);
		double x1=n(e,1,0);
		double x2=n(e,2,0);
		double x3=n(e,3,0);
		double y0=n(e,0,1);
		double y1=n(e,1,1);
		double y2=n(e,2,1);
		double y3=n(e,3,1);
		curl(e,0,0) = 0.5 * (-x1 + x3); curl(e,0,1) = 0.5 * (-y1 + y3);
		curl(e,1,0) = 0.5 * (x0 - x2);  curl(e,1,1) = 0.5 *  (y0 - y2);
		curl(e,2,0) = 0.5 * (x1 - x3);  curl(e,2,1) = 0.5 *  (y1 - y3);
		curl(e,3,0) = 0.5 * (-x0 + x2); curl(e,3,1) = 0.5 * (-y0 + y2);
	}

	return curlPhi;

}


IntrepidData::IntrepidData(size_t numCells, size_t cubDegree, const CellTopology& cellType)
: numCells(numCells), cubDegree(cubDegree),num_cubature_points(-1), cellType(cellType)
{

	// Integration rule
	Intrepid::DefaultCubatureFactory<double> cubFactory;
	myCub = cubFactory.create(cellType, cubDegree);
	num_cubature_points = myCub->getNumPoints();

	/*
	 * Scalar values at quadrature points
	 */
	nuGP                = RCP<FieldContainer>(new FieldContainer(numCells, num_cubature_points));
	sigmaGP             = RCP<FieldContainer>(new FieldContainer(numCells, num_cubature_points));
	weighted_measure    = RCP<FieldContainer>(new FieldContainer(numCells, num_cubature_points));
	normB               = RCP<FieldContainer>(new FieldContainer(numCells, num_cubature_points));

	/*
	 * basis at quadrature points
	 */
	basis_at_cub_points = RCP<FieldContainer>(new FieldContainer(numFields, num_cubature_points));
	grad_at_cub_points  = RCP<FieldContainer>(new FieldContainer(numFields, num_cubature_points, spaceDim));

	// Allocate and format arrays
	cub_points =   RCP<FieldContainer>(new FieldContainer           (num_cubature_points, spaceDim));
	cub_weights =  RCP<FieldContainer>(new FieldContainer           (num_cubature_points));
	jacobian =     RCP<FieldContainer>(new FieldContainer (numCells, num_cubature_points, spaceDim, spaceDim));
	jacobian_inv = RCP<FieldContainer>(new FieldContainer (numCells, num_cubature_points, spaceDim, spaceDim));
	jacobian_det = RCP<FieldContainer>(new FieldContainer (numCells, num_cubature_points));

	/*
	 * Material properties multiplied by quadrature weights
	 */
	weighted_measureNu    = RCP<FieldContainer>(new FieldContainer (numCells, num_cubature_points));
	weighted_measureSigma = RCP<FieldContainer>(new FieldContainer (numCells, num_cubature_points));

	/*
	 * Function values
	 */
	transformedVal         = RCP<FieldContainer>(new FieldContainer (numCells, numFields, num_cubature_points));
	weightedTransformedVal = RCP<FieldContainer>(new FieldContainer (numCells, numFields, num_cubature_points));

	/*
	 * Function gradients
	 */
	transformedGrad         = RCP<FieldContainer>(new FieldContainer(numCells, numFields, num_cubature_points, spaceDim));
	weightedTransformedGrad = RCP<FieldContainer>(new FieldContainer(numCells, numFields, num_cubature_points, spaceDim));

	/*
	 * gradA -- gradient of vector potential (C,P,D)
	 */
	gradA = RCP<FieldContainer>(new FieldContainer(numCells, num_cubature_points, spaceDim));

	/*
	 * Initialize geometry independent data
	 */

	// Select HGRAD basis -- note that this is hardwired for "QUADS" but
	//   it should be possible to get rid of this assumption
	Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, FieldContainer > quadBasis;

	// Compute basis functions and gradients at quadrature points
	myCub->getCubature(*cub_points, *cub_weights);
	quadBasis.getValues(*basis_at_cub_points, *cub_points, Intrepid::OPERATOR_VALUE);
	quadBasis.getValues(*grad_at_cub_points,  *cub_points, Intrepid::OPERATOR_GRAD);

}


void IntrepidData::computeNormB
(
        const FieldContainer& cell_nodes,
        const FieldContainer& A3Np1,
        double B0[2]
)
{
        /*
         * PRIOR TO CALLING THIS FUNCTION:
         * FOLLOWING LINE NEEDS TO BE IMPLEMENTED SOMEWHERE
         */
        // FIFE_2D_Quadrature_Intrepid_ElementErrorChecking(npe,spaceDim,elementBlock);

        // Compute Jacobian
        Intrepid::CellTools<double>::setJacobian(*jacobian, *cub_points, cell_nodes, cellType);
        Intrepid::CellTools<double>::setJacobianInv(*jacobian_inv, *jacobian);
        Intrepid::CellTools<double>::setJacobianDet(*jacobian_det, *jacobian);


        /*
         * Following function needs to be called somewhere and some how;
         */
//      FIFE_2D_Quadrature_Intrepid_Assert_detJacobian(jacobian_det,numCubPoints,elementBlock,ep);

        // Multiply jacobian_det with cubature weights
        Intrepid::FunctionSpaceTools::computeCellMeasure<double>(*weighted_measure, *jacobian_det, *cub_weights);
        Intrepid::FunctionSpaceTools::HGRADtransformGRAD<double>(*transformedGrad,*jacobian_inv,*grad_at_cub_points);
        Intrepid::FunctionSpaceTools::evaluate<double>(*gradA,A3Np1,*transformedGrad);

        /*
         * Compute normB = norm(gradA) + B0
         * JAM: ARE WE GOING TO STICK WITH B0?
         */
        {
                double B0X=B0[0];
                double B0Y=B0[1];
                double gradX, gradY;
                for(int cell=0;cell<numCells;cell++){
                        for(int gp=0;gp<num_cubature_points;gp++){
                                gradX = gradA->operator()(cell,gp,0) + B0X;
                                gradY = gradA->operator()(cell,gp,1) + B0Y;
                                normB->operator()(cell,gp) = sqrt(gradX*gradX + gradY*gradY);
                        }

                }
        }

}


void IntrepidData::finish_quadrature
(
        FieldContainer& stiffness_matrices,
        FieldContainer& mass_matrices
)
{

	// Scale weighted_measure with material properties and constants
	Intrepid::FunctionSpaceTools::scalarMultiplyDataData<double>(*weighted_measureSigma,*sigmaGP,*weighted_measure,false);
	Intrepid::FunctionSpaceTools::scalarMultiplyDataData<double>(*weighted_measureNu,*nuGP,*weighted_measure,false);


	/*
	 * Mass matrix pre-cursors
	 */
	Intrepid::FunctionSpaceTools::HGRADtransformVALUE<double>(*transformedVal,*basis_at_cub_points);
	Intrepid::FunctionSpaceTools::multiplyMeasure<double>(*weightedTransformedVal,*weighted_measureSigma,*transformedVal);

	/*
	 * Stiffness matrix precursors
	 */
	Intrepid::FunctionSpaceTools::multiplyMeasure<double>(*weightedTransformedGrad,*weighted_measureNu,*transformedGrad);

	/*
	 * Initialize stiffness and mass matrices
	 */
	for(std::size_t cell=0;cell<numCells;cell++)
		for(std::size_t n=0;n<numFields;n++)
			for(std::size_t d=0;d<numFields;d++){
				stiffness_matrices(cell,n,d)=0;
				mass_matrices(cell,n,d)=0;
			}

	// Integrate to compute mass and remainder of stiffness matrices (note "SUM_INTO=false" on stiffness)
	bool SUM_INTO=false;
	Intrepid::FunctionSpaceTools::integrate<double>(mass_matrices,*transformedVal,*weightedTransformedVal,Intrepid::COMP_CPP);
	Intrepid::FunctionSpaceTools::integrate<double>(stiffness_matrices,*transformedGrad,*weightedTransformedGrad,Intrepid::COMP_CPP,SUM_INTO);


}

#if 0
/**
 * @param A3Np1 -- vector potential at step n+1 (C,F)
 */
void CartesianQuadrature
(
        const CellTopology& cellType,
        const FieldContainer& cell_nodes,
        const FieldContainer& A3Np1,
        double B0[2],
        FieldContainer& stiffness_matrices,
        FieldContainer& mass_matrices
){

        int spaceDim = cellType.getDimension();
        int numCells = cell_nodes.dimension(0);


        /*
         * PRIOR TO CALLING THIS FUNCTION:
         * FOLLOWING LINE NEEDS TO BE IMPLEMENTED IN MHD
         */
        // FIFE_2D_Quadrature_Intrepid_ElementErrorChecking(npe,spaceDim,elementBlock);

        // Integration rule
        Intrepid::DefaultCubatureFactory<double> cubFactory;
        int cubDegree=2;
        Teuchos::RCP<Intrepid::Cubature<double> > myCub = cubFactory.create(cellType, cubDegree);
        int numCubPoints = myCub->getNumPoints();

        // Select HGRAD basis
        Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, Intrepid::FieldContainer<double> > quadBasis;
        int numFields = quadBasis.getCardinality();

        // Allocate and format arrays
        Intrepid::FieldContainer<double> cub_points(numCubPoints, spaceDim);
        Intrepid::FieldContainer<double> cub_weights(numCubPoints);

        Intrepid::FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
        Intrepid::FieldContainer<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
        Intrepid::FieldContainer<double> jacobian_det(numCells, numCubPoints);
        Intrepid::FieldContainer<double> weighted_measure(numCells, numCubPoints);

        /*
         * Allow material properties to vary by gauss point
         */
        Intrepid::FieldContainer<double> nuGP(numCells, numCubPoints);
        Intrepid::FieldContainer<double> sigmaGP(numCells, numCubPoints);
        Intrepid::FieldContainer<double> weighted_measureNu(numCells, numCubPoints);
        Intrepid::FieldContainer<double> weighted_measureSigma(numCells, numCubPoints);

        /*
         * Function values
         */
        Intrepid::FieldContainer<double> basis_at_cub_points(numFields, numCubPoints);
        Intrepid::FieldContainer<double> transformedVal(numCells, numFields, numCubPoints);
        Intrepid::FieldContainer<double> weightedTransformedVal(numCells, numFields, numCubPoints);

        /*
         * Function gradients
         */
        Intrepid::FieldContainer<double> grad_at_cub_points(numFields, numCubPoints, spaceDim);
        Intrepid::FieldContainer<double> transformedGrad(numCells,numFields, numCubPoints, spaceDim);
        Intrepid::FieldContainer<double> weightedTransformedGrad(numCells,numFields, numCubPoints, spaceDim);

        /*
         * gradA -- gradient of vector potential (C,P,D)
         */
        Intrepid::FieldContainer<double> gradA(numCells, numCubPoints, spaceDim);
        Intrepid::FieldContainer<double> normB(numCells, numCubPoints);
        /*
         * Mass and stiffness matrices
         */
//      Intrepid::FieldContainer<double> stiffness_matrices(numCells, numFields, numFields);
//      Intrepid::FieldContainer<double> mass_matrices(numCells, numFields, numFields);

        // Compute basis functions and gradients at gauss points
        myCub->getCubature(cub_points, cub_weights);
        quadBasis.getValues(basis_at_cub_points, cub_points, Intrepid::OPERATOR_VALUE);
        quadBasis.getValues(grad_at_cub_points, cub_points, Intrepid::OPERATOR_GRAD);

        // Compute Jacobian
        Intrepid::CellTools<double>::setJacobian(jacobian, cub_points, cell_nodes, cellType);
        Intrepid::CellTools<double>::setJacobianInv(jacobian_inv, jacobian);
        Intrepid::CellTools<double>::setJacobianDet(jacobian_det, jacobian);

        /*
         * Following function needs to be called somewhere and some how;
         */
//      FIFE_2D_Quadrature_Intrepid_Assert_detJacobian(jacobian_det,numCubPoints,elementBlock,ep);

        // Multiply jacobian_det with cubature weights
        Intrepid::FunctionSpaceTools::computeCellMeasure<double>(weighted_measure, jacobian_det, cub_weights);

        /*
         * Jacobian Specific Terms:
         * 0) gradA
         * 1) gradA dot gradPhi
         * 2) norm_gradA
         */
        Intrepid::FunctionSpaceTools::HGRADtransformGRAD<double>(transformedGrad,jacobian_inv,grad_at_cub_points);
        Intrepid::FunctionSpaceTools::evaluate<double>(gradA,A3Np1,transformedGrad);

        /*
         * Compute normB = norm(gradA) + B0
         * JAM: ARE WE GOING TO STICK WITH B0?
         */
        {
                double B0X=B0[0];
                double B0Y=B0[1];
                double gradX, gradY;
                for(int cell=0;cell<numCells;cell++){
                        for(int gp=0;gp<numCubPoints;gp++){
                                gradX = gradA(cell,gp,0) + B0X;
                                gradY = gradA(cell,gp,1) + B0Y;
                                normB(cell,gp) = sqrt(gradX*gradX + gradY*gradY);
                        }

                }
        }

        /*
         * CALL MagneticPropertyInterface
         * ** computes sigmaGP, numGP, and slopeNu
         */
//      matPropI.getProperties(normB,sigmaGP,nuGP);



        // Scale weighted_measure with 'massm' and 'diffc' for mass and stiffness matrices respectively
        Intrepid::FunctionSpaceTools::scalarMultiplyDataData<double>(weighted_measureSigma,sigmaGP,weighted_measure,false);
        Intrepid::FunctionSpaceTools::scalarMultiplyDataData<double>(weighted_measureNu,nuGP,weighted_measure,false);

        /*
         * Mass matrix pre-cursors
         */
        Intrepid::FunctionSpaceTools::HGRADtransformVALUE<double>(transformedVal,basis_at_cub_points);
        Intrepid::FunctionSpaceTools::multiplyMeasure<double>(weightedTransformedVal,weighted_measureSigma,transformedVal);

        /*
         * Stiffness matrix precursors
         */
        Intrepid::FunctionSpaceTools::multiplyMeasure<double>(weightedTransformedGrad,weighted_measureNu,transformedGrad);

        for(int cell=0;cell<numCells;cell++)
                for(int n=0;n<numFields;n++)
                        for(int d=0;d<numFields;d++){
                                stiffness_matrices(cell,n,d)=0;
                                mass_matrices(cell,n,d)=0;
                        }


        // Integrate to compute mass and remainder of stiffness matrices (note "SUM_INTO=false" on stiffness)
        bool SUM_INTO=false;
        Intrepid::FunctionSpaceTools::integrate<double>(mass_matrices,transformedVal,weightedTransformedVal,Intrepid::COMP_CPP);
        Intrepid::FunctionSpaceTools::integrate<double>(stiffness_matrices,transformedGrad,weightedTransformedGrad,Intrepid::COMP_CPP,SUM_INTO);

}

void CartesianJacobian
(
                const CellTopology& cellType,
                const FieldContainer& cell_nodes,
                const FieldContainer& A3Np1,
                double B0[2],
                FieldContainer& stiffness_matrices,
                FieldContainer& mass_matrices
){

        int spaceDim = cellType.getDimension();
        int numCells = cell_nodes.dimension(0);

        /*
         * PRIOR TO CALLING THIS FUNCTION:
         * FOLLOWING LINE NEEDS TO BE IMPLEMENTED IN MHD
         */
        // FIFE_2D_Quadrature_Intrepid_ElementErrorChecking(npe,spaceDim,elementBlock);

        // Integration rule
        Intrepid::DefaultCubatureFactory<double> cubFactory;
        int cubDegree=2;
        Teuchos::RCP<Intrepid::Cubature<double> > myCub = cubFactory.create(cellType, cubDegree);
        int numCubPoints = myCub->getNumPoints();

        // Select HGRAD basis
        Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, Intrepid::FieldContainer<double> > quadBasis;
        int numFields = quadBasis.getCardinality();

        // Allocate and format arrays
        Intrepid::FieldContainer<double> cub_points(numCubPoints, spaceDim);
        Intrepid::FieldContainer<double> cub_weights(numCubPoints);

        Intrepid::FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
        Intrepid::FieldContainer<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
        Intrepid::FieldContainer<double> jacobian_det(numCells, numCubPoints);
        Intrepid::FieldContainer<double> weighted_measure(numCells, numCubPoints);

        /*
         * Allow material properties to vary by gauss point
         */
        Intrepid::FieldContainer<double> nuGP(numCells, numCubPoints);
        Intrepid::FieldContainer<double> sigmaGP(numCells, numCubPoints);
        Intrepid::FieldContainer<double> weighted_measureNu(numCells, numCubPoints);
        Intrepid::FieldContainer<double> weighted_measureSigma(numCells, numCubPoints);

        /*
         * Function values
         */
        Intrepid::FieldContainer<double> basis_at_cub_points(numFields, numCubPoints);
        Intrepid::FieldContainer<double> transformedVal(numCells, numFields, numCubPoints);
        Intrepid::FieldContainer<double> weightedTransformedVal(numCells, numFields, numCubPoints);

        /*
         * Function gradients
         */
        Intrepid::FieldContainer<double> grad_at_cub_points(numFields, numCubPoints, spaceDim);
        Intrepid::FieldContainer<double> transformedGrad(numCells,numFields, numCubPoints, spaceDim);
        Intrepid::FieldContainer<double> weightedTransformedGrad(numCells,numFields, numCubPoints, spaceDim);

        /*
         * gradA -- gradient of vector potential (C,P,D)
         * gradA_dot_gradPhi -- gradient of vector potential dotted with gradient of shape functions (C,F,P)
         */
        Intrepid::FieldContainer<double> gradA(numCells, numCubPoints, spaceDim);
        Intrepid::FieldContainer<double> gradA_dot_gradPhi(numCells, numFields, numCubPoints);
        Intrepid::FieldContainer<double> weighted_gradA_dot_gradPhi(numCells, numFields, numCubPoints);
        Intrepid::FieldContainer<double> normB(numCells, numCubPoints);
        Intrepid::FieldContainer<double> slopeNu(numCells, numCubPoints);
        Intrepid::FieldContainer<double> slopeWt(numCells, numCubPoints);

        /*
         * Mass and stiffness matrices
         */
//      Intrepid::FieldContainer<double> stiffness_matrices(numCells, numFields, numFields);
//      Intrepid::FieldContainer<double> mass_matrices(numCells, numFields, numFields);

        // Compute basis functions and gradients at gauss points
        myCub->getCubature(cub_points, cub_weights);
        quadBasis.getValues(basis_at_cub_points, cub_points, Intrepid::OPERATOR_VALUE);
        quadBasis.getValues(grad_at_cub_points, cub_points, Intrepid::OPERATOR_GRAD);

        // Compute Jacobian
        Intrepid::CellTools<double>::setJacobian(jacobian, cub_points, cell_nodes, cellType);
        Intrepid::CellTools<double>::setJacobianInv(jacobian_inv, jacobian);
        Intrepid::CellTools<double>::setJacobianDet(jacobian_det, jacobian);

        /*
         * Following function needs to be called somewhere and some how;
         */
//      FIFE_2D_Quadrature_Intrepid_Assert_detJacobian(jacobian_det,numCubPoints,elementBlock,ep);

        // Multiply jacobian_det with cubature weights
        Intrepid::FunctionSpaceTools::computeCellMeasure<double>(weighted_measure, jacobian_det, cub_weights);

        /*
         * Jacobian Specific Terms:
         * 0) gradA
         * 1) gradA dot gradPhi
         * 2) norm_gradA
         */
        Intrepid::FunctionSpaceTools::HGRADtransformGRAD<double>(transformedGrad,jacobian_inv,grad_at_cub_points);
        Intrepid::FunctionSpaceTools::evaluate<double>(gradA,A3Np1,transformedGrad);
        Intrepid::FunctionSpaceTools::dotMultiplyDataField<double>(gradA_dot_gradPhi,transformedGrad,gradA);

        /*
         * Compute normB = norm(gradA) + B0
         */
        {
                double B0X=B0[0];
                double B0Y=B0[1];
                double gradX, gradY;
                for(int cell=0;cell<numCells;cell++){
                        for(int gp=0;gp<numCubPoints;gp++){
                                gradX = gradA(cell,gp,0) + B0X;
                                gradY = gradA(cell,gp,1) + B0Y;
                                normB(cell,gp) = sqrt(gradX*gradX + gradY*gradY);
                        }

                }
        }

        /*
         * CALL MagneticPropertyInterface
         * ** computes sigmaGP, numGP, and slopeNu
         */
//      matPropI.getProperties(normB,sigmaGP,nuGP,slopeNu);

        // Scale weighted_measure with 'massm' and 'diffc' for mass and stiffness matrices respectively
        Intrepid::FunctionSpaceTools::scalarMultiplyDataData<double>(weighted_measureSigma,sigmaGP,weighted_measure,false);
        Intrepid::FunctionSpaceTools::scalarMultiplyDataData<double>(weighted_measureNu,nuGP,weighted_measure,false);

        /*
         * Mass matrix pre-cursors
         */
        Intrepid::FunctionSpaceTools::HGRADtransformVALUE<double>(transformedVal,basis_at_cub_points);
        Intrepid::FunctionSpaceTools::multiplyMeasure<double>(weightedTransformedVal,weighted_measureSigma,transformedVal);

        /*
         * Stiffness matrix precursors
         */

        Intrepid::FunctionSpaceTools::multiplyMeasure<double>(weightedTransformedGrad,weighted_measureNu,transformedGrad);

        /*
         * FINAL Steps for Jacobian Calculation
         * 0) Given normB, compute slope of reluctivity curve: slopeNu
         * 1) Compute slopeWt = slopeNu * weighted_measure / normB
         * 2) compute weighted_gradA_dot_gradPhi
         * 3) integrate (weighted_gradA_dot_gradPhi,gradA_dot_gradPhi)
         */


        // 1)
        {
                for(int cell=0;cell<numCells;cell++){
                        for(int gp=0;gp<numCubPoints;gp++){
                                slopeWt(cell,gp) = weighted_measure(cell,gp) * slopeNu(cell,gp) / normB(cell,gp);
                        }

                }
        }
        // 2)
        Intrepid::FunctionSpaceTools::multiplyMeasure<double>(weighted_gradA_dot_gradPhi,slopeWt,gradA_dot_gradPhi);

        /**
         * 3)
         * 3a) Re-Initialize stiffness and mass matrices
         * 3b) Integrate
         * */

        // 3a)
        for(int cell=0;cell<numCells;cell++)
                for(int n=0;n<numFields;n++)
                        for(int d=0;d<numFields;d++){
                                stiffness_matrices(cell,n,d)=0;
                                mass_matrices(cell,n,d)=0;
                        }
        // 3b)
        Intrepid::FunctionSpaceTools::integrate<double>(stiffness_matrices,gradA_dot_gradPhi,weighted_gradA_dot_gradPhi,Intrepid::COMP_CPP);

        // Integrate to compute mass and remainder of stiffness matrices (note "SUM_INTO=true" on stiffness)
        bool SUM_INTO=true;
        Intrepid::FunctionSpaceTools::integrate<double>(mass_matrices,transformedVal,weightedTransformedVal,Intrepid::COMP_CPP);
        Intrepid::FunctionSpaceTools::integrate<double>(stiffness_matrices,transformedGrad,weightedTransformedGrad,Intrepid::COMP_CPP,SUM_INTO);

}


//void assertProperties(double conductivity, double reluctivity, int element_gid) const{
//
//      double sigma = conductivity;
//      double nu = reluctivity;
//
//      if ( (!finite(sigma) || sigma < 0.) ) {
//              std::ostringstream oss;
//              oss << "Mag::FIFE_2D_Quadrature_Intrepid(): "
//                              << "Scalar conductivity is < zero or trashed! (processor,block,global number,sigma)=("
//                              << rank << ","
//                              << block_id << ","
//                              << element_gid << ","
//                              << sigma << ")";
//              event_handler->log_error(oss.str());
//      }
//
//      if ( (!finite(nu) || nu <= 0. ) ) {
//              std::ostringstream oss;
//              oss <<"Mag::FIFE_2D_Quadrature_Intrepid(): "
//                              << "Scalar reluctivity is <= zero or trashed! (processor,block,global number,nu)=("
//                              << rank << ","
//                              << block_id << ","
//                              << element_gid << ","
//                              << nu << ")";
//              event_handler->log_error(oss.str());
//      }
//}

//      sigma = Scalar_Conductivity(ep);
//      nu    = Scalar_Reluctivity(ep);
//      FIFE_2D_Quadrature_Intrepid_AssertMaterialProperties(sigma,nu,elementBlock,ep);
//      if ( solve_mag_equation == A_EQUATION ) {
//              massm = kappa3*sigma;
//              diffc = nu/kappa5;
//      } else { // B_EQUATION
//              massm = 1.;
//              diffc = nu/(sigma*kappa3*kappa5);
//      }
//      for(int gp=0;gp<numCubPoints;gp++){
//              nuGP(cell,gp) = diffc;
//              sigmaGP(cell,gp) = massm;
//      }

//void Assert_detJacobian
//(
//              FieldContainer& jacobianDet,
//              int numCubPoints
//)
//{
//
//      /*
//       * THIS FUNCTION APPLIES TO '1' CELL / ELEMENT
//       */
//      int cell = 0;
//      for(int gp=0;gp<numCubPoints;gp++){
//              double detJacobian = jacobianDet(cell,gp);
//              if ( (!finite(detJacobian) || detJacobian <= 0.) ) {
//                      ostringstream oss;
//                      oss<<"Mag::FIFE_2D_Quadrature_Intrepid(): "
//                                      << "Element inverted or trashed! (processor,block,global number,gauss point,detJ)=("
//                                      << comm.rank() << ","
//                                      << block->Id() << ","
//                                      << element->globalNumber() << ","
//                                      << gp << ","
//                                      << detJacobian << ")";
//                      event_handler->log_error(oss.str());
//              }
//      }
//}
#endif

} // FIFE_2D_Quadrature


//#endif // HAVE_INTREPID
//
//#endif // TWO_D
//
//#else
//
//// some compilers will not accept an empty file for compilation
//// a nop is provided for this case
//
//void FIFE_2D_IntrepidQuadrature_noop() { return;}
//
//#endif // SNL_MHD
