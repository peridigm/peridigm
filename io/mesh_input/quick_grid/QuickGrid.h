/*
 * QuickGrid.h
 *
 *  Created on: Mar 12, 2011
 *      Author: jamitch
 */

#ifndef QUICKGRID_H_
#define QUICKGRID_H_

#include <vector>
#include "utilities/Array.h"
#include "utilities/Vector.h"
#include "QuickGridData.h"
#include <functional>


namespace QUICKGRID {

using std::size_t;
using std::tr1::shared_ptr;
using UTILITIES::Array;
using UTILITIES::Vector3D;
using std::binary_function;

bool NoOpNormFunction (const double* u, const double* v, double r);
bool SphericalNormFunction (const double* u, const double* v, double r);

typedef bool(*NormFunctionPointer)(const double*,const double*, double);
const NormFunctionPointer SphericalNorm = SphericalNormFunction;
const NormFunctionPointer NoOpNorm = NoOpNormFunction;

class Spec1D;
class SpecRing2D;
class QuickGridMeshGenerationIterator;

Array<double> getDiscretization(const Spec1D& spec);
Array<double> getDiscretization(const Spec1D& xSpec, const Spec1D& ySpec);
Array<double> getDiscretization(const Spec1D& xSpec, const Spec1D& ySpec, const Spec1D& zSpec);
Array<double> getDiscretization(const SpecRing2D& spec);
Array<double> getDiscretization(const SpecRing2D& spec, const Spec1D& axisSpec);
QuickGridData getDiscretization(size_t rank, QuickGridMeshGenerationIterator &cellIter);
QuickGridData allocatePdGridData(size_t numCells, size_t dimension);
shared_ptr<QuickGridMeshGenerationIterator> getMeshGenerator(size_t numProcs, const std::string& json_filename);
void print_meta_data(const QuickGridData& gridData, const std::string& label="");


class Horizon {
private:
	/*
	 * This is one-sided horizon
	 */
	size_t horizon;
	size_t nx;
public:
	explicit Horizon(size_t h, size_t numCells): horizon(h), nx(numCells) {}
	/**
	 * This function always includes the "ith" cell in its count
	 */
	size_t numCells(size_t i) const { return end(i)-start(i); }
	/*
	 * This is the first cell in the neighbor hood --
	 * Note that it could be the ith cell if "i" is on the boundary;
	 * In conjunction with 'numCells(i), this function can be used to
	 * generate the list of cells in the neighborhood; Like this:
	 * for(int j=start(i);j<start(i)+numCells(i);j++) {...}
	 * Note that this neighborhood will include "i"
	 */
	size_t start(size_t i) const { return i > horizon ? i-horizon : 0; }

private:
	size_t end(size_t i) const { return i<(nx-horizon) ? i+horizon+1 : nx; }

};

class HorizonIterator {
public:
	virtual ~HorizonIterator() {}
	virtual bool hasNextCell() const = 0;
 	virtual size_t numCells() const = 0;
	virtual size_t nextCell()  = 0;
};


class RingHorizon {
private:
	/*
	 * This is one-sided horizon
	 */
	size_t horizon;
	/*
	 * Number of cells in ring
	 */
	size_t numRays;

public:
	explicit RingHorizon(size_t h, size_t numCells): horizon(h), numRays(numCells) {}
	class RingHorizonIterator;
	friend class RingHorizonIterator;
	class RingHorizonIterator : public HorizonIterator {
	public:
		explicit RingHorizonIterator(size_t i,  size_t h, size_t nR) : HorizonIterator(), horizon(h), numRays(nR), ray(i), nextRay(start(i)), rayIter(0) {}
		virtual ~RingHorizonIterator(){}
		/*
		 * Is there another cell in the horizon
		 */
		bool hasNextCell() const { return rayIter<numCells(); }
		/*
		 * Number of cells in the horizon -- includes "this (ray)" cell
		 */
		size_t numCells() const { return 2*horizon+1; }
		/*
		 * Next cell in the horizon; returns -1 if function was improperly called
		 */
		size_t nextCell()  { int tmp=-1; if(hasNextCell()) tmp = nextRay; nextRay=incrementNextRay(); return tmp; }
	private:
		// Data
		size_t horizon, numRays;
		size_t ray, nextRay, rayIter;
		// Methods
		size_t start(size_t i) const { return i >= horizon ? i-horizon : numRays-(horizon-i); }
		size_t incrementNextRay()  { rayIter++; return (nextRay < numRays-1) ? ++nextRay : 0 ; }

	};

	/**
	 * This function always includes the "ith" cell in its count
	 */
	RingHorizonIterator horizonIterator(size_t i) const { return RingHorizonIterator(i, this->horizon,this->numRays); }

};

class Spec1D {
public:
	explicit Spec1D(size_t n, double x0, double L) : numCells(n), xStart(x0), xLength(L){}
	size_t getNumCells() const { return numCells; }
	double getX0() const { return xStart; }
	double getLength() const { return xLength; }
	double getCellSize() const { return xLength/numCells; }
	Horizon getCellHorizon(double horizon) const;
	/* This function should be called on a spec constructed in radians */
	RingHorizon getRingCellHorizon(double horizon, double ringRadius) const;

private:
	size_t numCells;
	double xStart, xLength;
	size_t getCellNeighborhoodSize(double horizon,double ringRadius=1.0) const;

};

/**
 *
 */
class SpecRing2D {
public:
	explicit SpecRing2D(Vector3D center, double innerRadius, double outerRadius, int numRings);
	Vector3D getCenter() const { return c; }
	double getrI() const { return rI; }
	double getr0() const { return r0; }
	double getRingThickness() const;
	size_t getNumRings() const { return numRings; }
	size_t getNumRays() const { return numRays; }
	size_t getNumCells() const { return numCells; }
	Spec1D getRingSpec() const {return Spec1D(getNumRays(), 0.0, 2.0*M_PI); }
	Spec1D getRaySpec() const  {return Spec1D(getNumRings(), getrI(), getRingThickness()); }


private:
	// This is just a 3d point so default copy is not costly unless
	//  we are creating thousands of SpecRing2D objects
	Vector3D c;
	double rI, r0;
	size_t numRings, numRays, numCells;
};


struct Cell3D {
public:
	Cell3D(size_t ii, size_t jj, size_t kk) : i(ii), j(jj), k(kk) {}
	size_t i, j, k;
};


class QuickGridMeshGenerationIterator {

public:
	QuickGridMeshGenerationIterator(size_t numProcs, NormFunctionPointer norm=NoOpNorm);
	virtual QuickGridData allocatePdGridData() const = 0;
	virtual std::pair<Cell3D,QuickGridData> beginIterateProcs(QuickGridData& pdGridDataProc0);
	bool hasNextProc() const { return iteratorProc < numProcs ? true : false; }
	size_t proc() const { return iteratorProc; }
	virtual std::pair<Cell3D,QuickGridData>  nextProc(Cell3D cellLocator, QuickGridData& pdGridDataProcN);
	virtual size_t getNumGlobalCells() const = 0;
	virtual size_t getDimension() const = 0;
	virtual ~QuickGridMeshGenerationIterator() {}
	virtual std::pair<Cell3D,QuickGridData> computePdGridData(size_t proc, Cell3D cellLocator, QuickGridData& pdGridData, NormFunctionPointer norm = NoOpNorm) const = 0;
	static std::vector<size_t> getNumCellsPerProcessor(size_t globalNumCells, size_t numProcs);
	int getNumProcs() const { return numProcs; }

private:
	int iteratorProc, numProcs;
	NormFunctionPointer neighborHoodNorm;
};

class TensorProduct3DMeshGenerator : public QuickGridMeshGenerationIterator {

public:
	explicit TensorProduct3DMeshGenerator
	(
			size_t numProc,
			double horizonRadius,
			const Spec1D& xSpec,
			const Spec1D& ySpec,
			const Spec1D& zSpec,
			NormFunctionPointer norm=NoOpNorm
	);
	virtual ~TensorProduct3DMeshGenerator() {}
	virtual QuickGridData allocatePdGridData() const;
	size_t getNumGlobalCells() const { return globalNumberOfCells; }
	size_t getDimension() const { return 3; }
	std::pair<Cell3D,QuickGridData> computePdGridData(size_t proc, Cell3D cellLocator, QuickGridData& pdGridData, NormFunctionPointer norm = NoOpNorm) const;
private:
	double horizonRadius;
	size_t globalNumberOfCells;
	std::vector<Spec1D> specs;
	std::vector<Horizon> H;
	std::vector<size_t> cellsPerProc;
	// This is for computing the tensor product space
	Array<double> xx;
	Array<double> yy;
	Array<double> zz;

	size_t getSizeNeighborList(size_t proc, Cell3D cellLocator) const;
	size_t computeNumNeighbors(size_t i, size_t j, size_t k) const;

};

class TensorProductCylinderMeshGenerator : public QuickGridMeshGenerationIterator {

public:
	explicit TensorProductCylinderMeshGenerator
	(
			size_t numProcs,
			double radiusHorizon,
			const SpecRing2D& rSpec, const
			Spec1D& zSpec,
			NormFunctionPointer norm=NoOpNorm
	);
	virtual ~TensorProductCylinderMeshGenerator() {}
	virtual QuickGridData allocatePdGridData() const;
	size_t getNumGlobalCells() const { return globalNumberOfCells; }
	size_t getDimension() const { return 3; }
	std::pair<Cell3D,QuickGridData> computePdGridData(size_t proc, Cell3D cellLocator, QuickGridData& pdGridData, NormFunctionPointer norm = NoOpNorm) const;
	const std::vector<Spec1D>& getTensorProductSpecs() const { return specs; }
	/*
	 * These are only public for testing purposes -- don't call these functions
	 */
	size_t getSizeNeighborList(size_t proc, Cell3D cellLocator) const;
	size_t computeNumNeighbors(size_t i, size_t j, size_t k) const;

private:
	// Data
	double horizonRadius;
	size_t globalNumberOfCells;
	SpecRing2D ringSpec;
	std::vector<Spec1D> specs;
	RingHorizon ringHorizon;
	std::vector<Horizon> H;
	std::vector<size_t> cellsPerProc;

	// This is for computing the tensor product space
	Array<double> rPtr;
	Array<double> thetaPtr;
	Array<double> zPtr;


};

class TensorProductSolidCylinder : public QuickGridMeshGenerationIterator {

public:
	explicit TensorProductSolidCylinder
	(
			size_t nProcs,
			double radius,
			size_t numRings,
			double zStart,
			double cylinderLength
	);
	virtual ~TensorProductSolidCylinder() {}
	virtual QuickGridData allocatePdGridData() const;
	size_t getNumGlobalCells() const { return globalNumberOfCells; }
	size_t getDimension() const { return 3; }
	std::pair<Cell3D,QuickGridData> computePdGridData(size_t proc, Cell3D cellLocator, QuickGridData& pdGridData, NormFunctionPointer norm = NoOpNorm) const;
	const std::vector<Spec1D>& getTensorProductSpecs() const { return specs; }

private:
	// Data
	double coreRadius, outerRadius, coreCellVolume;
	size_t globalNumberOfCells;
	std::vector<Spec1D> specs;
	std::vector<size_t> cellsPerProc;

	// This is for computing the tensor product space
	Array<double> rPtr;
	Array<double> thetaPtr;
	Array<double> zPtr;


};



}

#endif /* QUICKGRID_H_ */
