/*
 * quadrulator.h
 *
 *  Created on: Nov 13, 2010
 *      Author: awesome
 */

#ifndef QUADRULATOR_H_
#define QUADRULATOR_H_

#include <cstdio>
#include <Teuchos_ArrayRCP.hpp>

using std::size_t;
using Teuchos::ArrayRCP;

class Quadrulator {
public:
	Quadrulator(double x0, double y0, std::size_t nx, double hx, std::size_t ny, double hy);
	ArrayRCP<size_t> getVertexLinks() { return vLinks; }
	ArrayRCP<double> getCoordinates() { return coord; }
	ArrayRCP<double> getElementCoordinates() { return localCoord; }
	size_t getNumVertices() const { return numVertices; }
	size_t getNpe() const { return npe; }
	size_t getNem() const { return numQuad; }
private:
	void computeVertexCoordinates();
	void computeLocalVertexCoordinates();
	void computeVertexLinks();

private:
	size_t nx, ny, npe, numVertices, numQuad;
	double x0, y0, hx, hy;
	ArrayRCP<double> coord, localCoord;
	ArrayRCP<size_t> vLinks;
};



#endif /* QUADRULATOR_H_ */
