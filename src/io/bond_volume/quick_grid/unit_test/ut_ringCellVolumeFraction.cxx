//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include <string>
#include "QuickGrid.h"
#include "QuickGridData.h"
#include "calculators.h"
#include "NeighborhoodList.h"
#include "PdZoltan.h"
#include "PdutMpiFixture.h"
#include "Field.h"
#include "Vector3D.h"
#include "Array.h"
#include <set>
#include <Teuchos_UnitTestHarness.hpp>
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include <Teuchos_ParameterList.hpp>
#ifdef USE_YAML
  #include <Teuchos_YamlParameterListCoreHelpers.hpp>
#endif

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <iostream>
#include <cmath>
using namespace Pdut;

using namespace QUICKGRID;
using UTILITIES::Vector3D;
using UTILITIES::Array;
using Field_NS::Field;
using Field_NS::FieldSpec;
using std::shared_ptr;

using std::size_t;
using std::string;
using std::cout;

static size_t myRank;
static size_t numProcs;
const string yaml_file_name="./input_files/ut_ringCellVolumeFraction.yaml";

void compute_neighborhood_volumes
(
		const PDNEIGH::NeighborhoodList& list,
		Field<double>& neighborhoodVol,
		Field<double>& naiveNeighborhoodVol,
		Array<double>& overlapCellVol,
		shared_ptr<double> xOverlapPtr,
		const BOND_VOLUME::QUICKGRID::RingVolumeFractionCalculator& calculator
)
{
	size_t N = list.get_num_owned_points();
	

	const int *neighPtr = list.get_local_neighborhood().get();
	const double *xOwned = list.get_owned_x().get();
	const double *xOverlap = xOverlapPtr.get();
	const double *cellVolOverlap = overlapCellVol.get();
	double *neighVol = neighborhoodVol.get();
	double *naiveNeighVol = naiveNeighborhoodVol.get();
	for(size_t p=0;p<N;p++, xOwned +=3, neighVol++,naiveNeighVol++){
		int numNeigh = *neighPtr; neighPtr++;
		/*
		 * initialize neighborhood to zero;
		 * computes contributions from neighborhood and does not
		 * include self volume
		 */
		*neighVol = 0.0;
		*naiveNeighVol = 0.0;

		const double *P = xOwned;
		/*
		 * Loop over neighborhood of point P and compute
		 * fractional volume
		 */
		for(int n=0;n<numNeigh;n++,neighPtr++){
			int localId = *neighPtr;
			const double *Q = &xOverlap[3*localId];
			double cellVolume = cellVolOverlap[localId];
			*neighVol += calculator(P,Q);
			*naiveNeighVol += cellVolume;
		}
	}
}

void compute_cell_volumes
(
		const PDNEIGH::NeighborhoodList& list,
		Field<double>& specialCellVolume,
		shared_ptr<double> xOverlapPtr,
		const BOND_VOLUME::QUICKGRID::RingVolumeFractionCalculator& calculator
)
{
	size_t N = list.get_num_owned_points();
	
	const double *xOwned = list.get_owned_x().get();
	double *vOwned = specialCellVolume.get();
	for(size_t p=0;p<N;p++, xOwned +=3, vOwned++){

		const double *P = xOwned;
		/*
		 * compute cell volume using quick grid quadrature
		 */
		*vOwned=calculator.cellVolume(P);

	}
}

TEUCHOS_UNIT_TEST(RingCellVolumeFraction, SphereTest) {

  Teuchos::ParameterList params;
  Teuchos::Ptr<Teuchos::ParameterList> params_ptr(&params);
  if (!QUICKGRID::file_exists(yaml_file_name)) {
    std::string msg = "\n\n**** Error, failed to open file " + yaml_file_name;
    msg += " (unit tests may need to be executed from within the directory where they reside).";
    std::cout << msg << std::endl;;
    throw std::runtime_error(msg);
  }
#ifdef USE_YAML
  Teuchos::updateParametersFromYamlFile(yaml_file_name, params_ptr);
#else
  std::string msg = "\n**** Error, unit test requires YAML.";
  std::cout << "\n" << msg << std::endl;;
  throw std::runtime_error(msg);
#endif

  Teuchos::ParameterList disc_params = params.sublist("Discretization");
  double horizon = disc_params.get<double>("Horizon");
  std::string type = disc_params.get<std::string>("Type");

  /*
	 * Construct ring spec
	 */
	// set ring at center point
	// Note that ring is only translated to the location
  Vector3D center;
  Teuchos::ParameterList mesh_params = params.sublist("QuickGrid").sublist("TensorProductCylinderMeshGenerator");
  center[0] = mesh_params.get<double>("Ring Center X");
  center[1] = mesh_params.get<double>("Ring Center Y");
  center[2] = mesh_params.get<double>("Ring Center Z");

  double innerRadius = mesh_params.get<double>("Inner Radius");
	TEST_ASSERT(16.0==innerRadius);
	double outerRadius = mesh_params.get<double>("Outer Radius");
	TEST_ASSERT(17.0==outerRadius);
	double cylinderLength=mesh_params.get<double>("Cylinder Length");
	TEST_ASSERT(1.0==cylinderLength);
	/*
	 * Note that zStart is set to the z - coordinate of 'center'
	 */
	double zStart=center[2];
	TEST_ASSERT(0.0==zStart);
	size_t numRings = mesh_params.get<int>("Number Points Radius");
	SpecRing2D ring2dSpec(center,innerRadius,outerRadius,numRings);

	// Create 1d Spec along cylinder axis
	// Compute number of cells along length of cylinder so that aspect ratio
	// is cells is approximately 1.
	// Cell sizes along axis are not exactly "cellSize" since last cell
	// would be a fraction of a cellSize -- so 1 is added to numCellsAlongAxis.
	// Actual cell sizes are slightly smaller than "cellSize" because of this.
	double cellSize = ring2dSpec.getRaySpec().getCellSize();
	int numCellsAxis = (int)(cylinderLength/cellSize)+1;
	Spec1D axisSpec(numCellsAxis,zStart,cylinderLength);

	BOND_VOLUME::QUICKGRID::RingVolumeFractionCalculator calculator(ring2dSpec,axisSpec,horizon);

	// Create decomposition iterator
	QUICKGRID::TensorProductCylinderMeshGenerator cellPerProcIter(numProcs, horizon,ring2dSpec, axisSpec);
	QUICKGRID::QuickGridData gridData = QUICKGRID::getDiscretization(myRank, cellPerProcIter);;


	// This load-balances
	gridData = PDNEIGH::getLoadBalancedDiscretization(gridData);

	/*
	 * Create neighborhood with an enlarged horizon
	 */
	double cell_diagonal=calculator.get_cell_diagonal();
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon+cell_diagonal);
	Array<double> xOverlapArray;
	Array<double> vOverlapArray;
	{
		/*
		 * mesh coordinates overlap
		 */
		Epetra_BlockMap ownedMap(*list.getOwnedMap(3));
		Epetra_BlockMap overlapMap(*list.getOverlapMap(3));
		xOverlapArray = Array<double>(overlapMap.NumMyElements()*3);
		Epetra_Import importNDF(overlapMap,ownedMap);
		Epetra_Vector xOverlap(View,overlapMap,xOverlapArray.get());
		Epetra_Vector xOwned(View,ownedMap,list.get_owned_x().get());
		xOverlap.Import(xOwned,importNDF,Insert);
	}

	{
		/*
		 * volume overlap
		 */
		Epetra_BlockMap ownedMap(*list.getOwnedMap(1));
		Epetra_BlockMap overlapMap(*list.getOverlapMap(1));
		vOverlapArray = Array<double>(overlapMap.NumMyElements()*1);
		vOverlapArray.set(0.0);
		Epetra_Import importNDF(overlapMap,ownedMap);
		Epetra_Vector vOverlap(View,overlapMap,vOverlapArray.get());
		Epetra_Vector vOwned(View,ownedMap,gridData.cellVolume.get());
		vOverlap.Import(vOwned,importNDF,Insert);
	}





	/*
	 * Compute volume on neighborhood for every point in mesh;
	 * Points with a spherical neighborhood that are completely enclosed
	 * should have a volume that very closely matches the analytical value for a sphere
	 */
	FieldSpec neighVolSpec
	(Field_ENUM::TYPE_UNDEFINED,Field_ENUM::ELEMENT,Field_ENUM::SCALAR, Field_ENUM::CONSTANT,"neighVol");

	FieldSpec naiveNeighVolSpec
	(Field_ENUM::TYPE_UNDEFINED,Field_ENUM::ELEMENT,Field_ENUM::SCALAR,Field_ENUM::CONSTANT,"naiveNeighVol");

	FieldSpec quadratureCellVolSpec
	(Field_ENUM::TYPE_UNDEFINED,Field_ENUM::ELEMENT,Field_ENUM::SCALAR,Field_ENUM::CONSTANT,"quadratureCellVol");

	Field<double> neighVol(neighVolSpec,list.get_num_owned_points());
	Field<double> naiveNeighVol(naiveNeighVolSpec,list.get_num_owned_points());
	Field<double> quadratureCellVol(quadratureCellVolSpec,list.get_num_owned_points());
	Field<double> cellVol(Field_NS::VOLUME,gridData.cellVolume,list.get_num_owned_points());

        size_t N = list.get_num_owned_points();
	TEST_ASSERT(quadratureCellVol.get_num_points()==N);
	compute_cell_volumes(list,quadratureCellVol,xOverlapArray.get_shared_ptr(),calculator);
        
        TEST_ASSERT(neighVol.get_num_points()==N);
	TEST_ASSERT(naiveNeighVol.get_num_points()==N);
	compute_neighborhood_volumes(list,neighVol,naiveNeighVol,vOverlapArray,xOverlapArray.get_shared_ptr(),calculator);
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

	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}

