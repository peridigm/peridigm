#ifndef FIFE_2D_INTREPIDQUADRATURE_H_
#define FIFE_2D_INTREPIDQUADRATURE_H_

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_Cubature.hpp>
#include <Shards_CellTopology.hpp>
#include <Teuchos_RCP.hpp>

using Teuchos::ArrayRCP;
using Teuchos::RCP;


namespace FIFE_2D_IntrepidQuadrature {

typedef Intrepid::FieldContainer<double> FieldContainer;

typedef shards::CellTopology CellTopology;

typedef Intrepid::Cubature<double> Cubature;

/*
 * @param cell_nodes
 * Dimension: C,F,D
 *     numCells, numFields, spaceDim
 * @ param cellType
 * Currently only supports Quadrilateral -- perhaps we can template this parameter
 * @return cell area
 * Dimension: C
 *     numCells
 */
RCP<FieldContainer> computeCellArea(const FieldContainer& cell_nodes, const CellTopology& cellType);

/*
 * @param cell_nodes
 * Dimension: C,F,D
 *     numCells,numFields,spacDim
 * @ param cellType
 * Currently only supports Quadrilateral -- perhaps we can template this parameter
 * @return cell area
 * Dimension: C, F, D
 *     nnumCells, numFields, spaceDim
 */
RCP<FieldContainer> meanQuadratureCurlPhi(const FieldContainer& cell_nodes, const CellTopology& cellType);

class IntrepidData {

public:

        explicit IntrepidData(std::size_t numCells, std::size_t cubDegree, const CellTopology& cellType);
       // ~IntrepidData() {}
        std::size_t getNumCells() const { return numCells; }
        std::size_t getNumCubaturePoints() const { return num_cubature_points; }
        FieldContainer& getReluctivity() { return *nuGP; }
        FieldContainer& getConductivity()  { return *sigmaGP; }
        FieldContainer& getWeightedMeasure() { return *weighted_measure; }
        FieldContainer& getBasis() { return *basis_at_cub_points; }
        FieldContainer& getTransformedGrad() { return *transformedGrad; }
        FieldContainer& getNormB() { return *normB; }
        void computeNormB(
                        const FieldContainer& cell_nodes,
                        const FieldContainer& A3Np1,
                        double B0[2]
        );
        void finish_quadrature
        (
                FieldContainer& stiffness_matrices,
                FieldContainer& mass_matrices
        );


private:
        std::size_t numCells, cubDegree, num_cubature_points;

        const shards::CellTopology cellType;

        Teuchos::RCP<Cubature> myCub;
        /*
         * C,P
         * numCells, num_cubature_points
         */
        RCP<FieldContainer> nuGP, sigmaGP, weighted_measure, normB;
        /*
         * F,P
         * numFields, num_cubature_points
         */
        RCP<FieldContainer> basis_at_cub_points;

        /*
         * C,F,P,D
         * numCells, numFields, num_cubature_points, spaceDim
         */
        RCP<FieldContainer> transformedGrad;

        /*
         * P,D
         */
        RCP<FieldContainer> cub_points;

        /*
         * P
         */
        RCP<FieldContainer> cub_weights;

        /*
         * C,P,D,D
         */
        RCP<FieldContainer> jacobian, jacobian_inv;

        /*
         * C,P
         * numCells, num_cubature_points
         */
        RCP<FieldContainer> jacobian_det, weighted_measureNu, weighted_measureSigma;

        /*
         * C, F, P
         */
        RCP<FieldContainer> transformedVal, weightedTransformedVal;

        /*
         * F, P, D
         */
        RCP<FieldContainer> grad_at_cub_points;

        /*
         * C, F, P, D
         */
        RCP<FieldContainer> weightedTransformedGrad;

        /*
         * gradA -- gradient of vector potential (C,P,D)
         */
        RCP<FieldContainer> gradA;

};



#if 0
void CartesianQuadrature
(
                const CellTopology& cellType,
                const FieldContainer& cell_nodes,
                const FieldContainer& A3Np1,
                double B0[2],
                FieldContainer& stiffness_matrices,
                FieldContainer& mass_matrices
);

void CartesianJacobian
(
                const CellTopology& cellType,
                const FieldContainer& cell_nodes,
                const FieldContainer& A3Np1,
                double B0[2],
                FieldContainer& stiffness_matrices,
                FieldContainer& mass_matrices
);

#endif
}

#endif /* FIFE_2D_INTREPIDQUADRATURE_H_ */
