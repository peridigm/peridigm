/*! \file PdQuickGrid.h */

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

#ifndef PDQUICKGRID_H_
#define PDQUICKGRID_H_
#include <vector>
#include<tr1/memory>
#include <valarray>
#include "PdGridData.h"
#include "Epetra_BlockMap.h"

class Epetra_Comm;

namespace PdQuickGrid {
using std::valarray;
/**
 * Utilities
 */
template<class T> struct Deleter{
	void operator()(T* d) {
		delete [] d;
	}
};



bool NoOpNormFunction (const double* u, const double* v, double r);
bool SphericalNormFunction (const double* u, const double* v, double r);

typedef bool(*NormFunctionPointer)(const double*,const double*, double);
const NormFunctionPointer SphericalNorm = SphericalNormFunction;
const NormFunctionPointer NoOpNorm = NoOpNormFunction;

PdGridData allocatePdGridData(int numCells, int dimension);

class PdQRing2d;
class PdQPointSet1d;

std::tr1::shared_ptr<double> getDiscretization(const PdQPointSet1d& spec);
std::tr1::shared_ptr<double> getDiscretization(const PdQPointSet1d& xSpec, const PdQPointSet1d& ySpec);
std::tr1::shared_ptr<double> getDiscretization(const PdQPointSet1d& xSpec, const PdQPointSet1d& ySpec, const PdQPointSet1d& zSpec);
std::tr1::shared_ptr<double> getDiscretization(const PdQRing2d& spec);
std::tr1::shared_ptr<double> getDiscretization(const PdQRing2d& spec, const PdQPointSet1d& axisSpec);
const Epetra_BlockMap getOverlapMap(const Epetra_Comm& comm, const PdGridData& gridData, int ndf);
const Epetra_BlockMap getOwnedMap(const Epetra_Comm& comm,const PdGridData& gridData, int ndf);
std::tr1::shared_ptr<int> getLocalOwnedIds(const PdGridData& gridData, const Epetra_BlockMap& overlapMap);
std::tr1::shared_ptr<int> getLocalNeighborList(const PdGridData& gridData, const Epetra_BlockMap& overlapMap);

class Horizon {
private:
	/*
	 * This is one-sided horizon
	 */
	int horizon;
	int nx;
public:
	explicit Horizon(int h, int numCells): horizon(h), nx(numCells) {}
	/**
	 * This function always includes the "ith" cell in its count
	 */
	int numCells(int i) const { return end(i)-start(i); }
	/*
	 * This is the first cell in the neighbor hood --
	 * Note that it could be the ith cell if "i" is on the boundary;
	 * In conjunction with 'numCells(i), this function can be used to
	 * generate the list of cells in the neighbor; Like this:
	 * for(int j=start(i);j<start(i)+numCells(i);j++) {...}
	 * Note that this neighborhood will include "i"
	 */
	int start(int i) const { return i > horizon ? i-horizon : 0; }

private:
	int end(int i) const { return i<(nx-horizon) ? i+horizon+1 : nx; }

};

class HorizonIterator {
public:
	virtual ~HorizonIterator() {}
	virtual bool hasNextCell() const = 0;
 	virtual int numCells() const = 0;
	virtual int nextCell()  = 0;
};


class RingHorizon {
private:
	/*
	 * This is one-sided horizon
	 */
	int horizon;
	/*
	 * Number of cells in ring
	 */
	int numRays;

public:
	explicit RingHorizon(int h, int numCells): horizon(h), numRays(numCells) {}
	class RingHorizonIterator;
	friend class RingHorizonIterator;
	class RingHorizonIterator : public HorizonIterator {
	public:
		explicit RingHorizonIterator(int i,  int h, int nR) : HorizonIterator(), horizon(h), numRays(nR), ray(i), nextRay(start(i)), rayIter(0) {}
		virtual ~RingHorizonIterator(){}
		/*
		 * Is there another cell in the horizon
		 */
		bool hasNextCell() const { return rayIter<numCells(); }
		/*
		 * Number of cells in the horizon -- includes "this (ray)" cell
		 */
		int numCells() const { return 2*horizon+1; }
		/*
		 * Next cell in the horizon; returns -1 if function was improperly called
		 */
		int nextCell()  { int tmp=-1; if(hasNextCell()) tmp = nextRay; nextRay=incrementNextRay(); return tmp; }
	private:
		// Data
		int horizon, numRays;
		int ray, nextRay, rayIter;
		// Methods
		int start(int i) const { return i >= horizon ? i-horizon : numRays-(horizon-i); }
		int incrementNextRay()  { rayIter++; return (nextRay < numRays-1) ? ++nextRay : 0 ; }

	};

	/**
	 * This function always includes the "ith" cell in its count
	 */
	RingHorizonIterator horizonIterator(int i) const { return RingHorizonIterator(i, this->horizon,this->numRays); }

};

class PdQPointSet1d{
public:
	explicit PdQPointSet1d(int n, double x0, double L) : numCells(n), xStart(x0), xLength(L){}
	int getNumCells() const { return numCells; }
	double getX0() const { return xStart; }
	double getLength() const { return xLength; }
	double getCellSize() const { return xLength/numCells; }
	Horizon getCellHorizon(double horizon) const;
	/* This function should be called on a spec constructed in radians */
	RingHorizon getRingCellHorizon(double horizon, double ringRadius) const;

private:
	int numCells;
	double xStart, xLength;
	int getCellNeighborhoodSize(double horizon,double ringRadius=1.0) const;

};

/**
 *
 */
class PdQRing2d {
public:
	explicit PdQRing2d(valarray<double> center, double innerRadius, double outerRadius, int numRings);
	valarray<double> getCenter() const { return c; }
	double getrI() const { return rI; }
	double getr0() const { return r0; }
	double getRingThickness() const;
	int getNumRings() const { return numRings; }
	int getNumRays() const { return numRays; }
	int getNumCells() const { return numCells; }
	PdQPointSet1d getRingSpec() const {return PdQPointSet1d(getNumRays(), 0.0, 2.0*M_PI); }
	PdQPointSet1d getRaySpec() const  {return PdQPointSet1d(getNumRings(), getrI(), getRingThickness()); }


private:
	// This is just a 3d point so default copy is not costly unless
	//  we are creating thousands of PdQRing2d objects
	valarray<double> c;
	double rI, r0;
	int numRings, numRays, numCells;
};


struct Cell3D {
public:
	Cell3D(int ii, int jj, int kk) : i(ii), j(jj), k(kk) {}
	int i, j, k;
};


class PdQuickGridMeshGenerationIterator {

public:
	PdQuickGridMeshGenerationIterator(int numProcs, NormFunctionPointer norm=NoOpNorm);
	virtual PdGridData allocatePdGridData() const = 0;
	virtual std::pair<Cell3D,PdGridData> beginIterateProcs(PdGridData& pdGridDataProc0);
	bool hasNextProc() const { return iteratorProc < numProcs ? true : false; }
	int proc() const { return iteratorProc; }
	virtual std::pair<Cell3D,PdGridData>  nextProc(Cell3D cellLocator, PdGridData& pdGridDataProcN);
	virtual int getNumGlobalCells() const = 0;
	virtual int getDimension() const = 0;
	virtual ~PdQuickGridMeshGenerationIterator() {}
	virtual std::pair<Cell3D,PdGridData> computePdGridData(int proc, Cell3D cellLocator, PdGridData& pdGridData, NormFunctionPointer norm = NoOpNorm) const = 0;
	static std::vector<int> getNumCellsPerProcessor(int globalNumCells, int numProcs);
	int getNumProcs() const { return numProcs; }

private:
	int iteratorProc, numProcs;
	NormFunctionPointer neighborHoodNorm;
};

class TensorProduct3DMeshGenerator : public PdQuickGridMeshGenerationIterator {

public:
	explicit TensorProduct3DMeshGenerator
	(
			int numProc,
			double horizonRadius,
			const PdQPointSet1d& xSpec,
			const PdQPointSet1d& ySpec,
			const PdQPointSet1d& zSpec,
			NormFunctionPointer norm=NoOpNorm
	);
	~TensorProduct3DMeshGenerator() {}
	PdGridData allocatePdGridData() const;
	int getNumGlobalCells() const { return globalNumberOfCells; }
	int getDimension() const { return 3; }
	std::pair<Cell3D,PdGridData> computePdGridData(int proc, Cell3D cellLocator, PdGridData& pdGridData, NormFunctionPointer norm = NoOpNorm) const;
private:
	double horizonRadius;
	int globalNumberOfCells;
	std::vector<PdQPointSet1d> specs;
	std::vector<Horizon> H;
	std::vector<int> cellsPerProc;
	// This is for computing the tensor product space
	std::tr1::shared_ptr<double> xx;
	std::tr1::shared_ptr<double> yy;
	std::tr1::shared_ptr<double> zz;

	int getSizeNeighborList(int proc, Cell3D cellLocator) const;
	int computeNumNeighbors(int i, int j, int k) const;

};

class TensorProductCylinderMeshGenerator : public PdQuickGridMeshGenerationIterator {

public:
	explicit TensorProductCylinderMeshGenerator
	(
			int numProcs,
			double radiusHorizon,
			const PdQRing2d& rSpec, const
			PdQPointSet1d& zSpec,
			NormFunctionPointer norm=NoOpNorm
	);
	~TensorProductCylinderMeshGenerator() {}
	PdGridData allocatePdGridData() const;
	int getNumGlobalCells() const { return globalNumberOfCells; }
	int getDimension() const { return 3; }
	std::pair<Cell3D,PdGridData> computePdGridData(int proc, Cell3D cellLocator, PdGridData& pdGridData, NormFunctionPointer norm = NoOpNorm) const;
	const std::vector<PdQPointSet1d>& getTensorProductSpecs() const { return specs; }
	/*
	 * These are only public for testing purposes -- don't call these functions
	 */
	int getSizeNeighborList(int proc, Cell3D cellLocator) const;
	int computeNumNeighbors(int i, int j, int k) const;

private:
	// Data
	double horizonRadius;
	int globalNumberOfCells;
	PdQRing2d ringSpec;
	std::vector<PdQPointSet1d> specs;
	RingHorizon ringHorizon;
	std::vector<Horizon> H;
	std::vector<int> cellsPerProc;

	// This is for computing the tensor product space
	std::tr1::shared_ptr<double> rPtr;
	std::tr1::shared_ptr<double> thetaPtr;
	std::tr1::shared_ptr<double> zPtr;


};

class TensorProductSolidCylinder : public PdQuickGridMeshGenerationIterator {

public:
	explicit TensorProductSolidCylinder
	(
			int nProcs,
			double radius,
			int numRings,
			double zStart,
			double cylinderLength
	);
	~TensorProductSolidCylinder() {}
	PdGridData allocatePdGridData() const;
	int getNumGlobalCells() const { return globalNumberOfCells; }
	int getDimension() const { return 3; }
	std::pair<Cell3D,PdGridData> computePdGridData(int proc, Cell3D cellLocator, PdGridData& pdGridData, NormFunctionPointer norm = NoOpNorm) const;
	const std::vector<PdQPointSet1d>& getTensorProductSpecs() const { return specs; }

private:
	// Data
	double coreRadius, outerRadius, coreCellVolume;
	int globalNumberOfCells;
	std::vector<PdQPointSet1d> specs;
	std::vector<int> cellsPerProc;

	// This is for computing the tensor product space
	std::tr1::shared_ptr<double> rPtr;
	std::tr1::shared_ptr<double> thetaPtr;
	std::tr1::shared_ptr<double> zPtr;


};

}
#endif /* PDQUICKGRID_H_ */
