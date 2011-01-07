/*
 * quadrulator.C
 *
 *  Created on: Nov 13, 2010
 *      Author: awesome
 */

#include "quadrulator.h"
#include <iostream>

Quadrulator::Quadrulator(double x0, double y0, std::size_t nx, double hx, std::size_t ny, double hy)
:  nx(nx), ny(ny), npe(4), numVertices((nx+1)*(ny+1)), numQuad(nx*ny), x0(x0), y0(y0), hx(hx), hy(hy),
   coord(3*(nx+1)*(ny+1),0.0), localCoord(3*4*nx*ny,0.0), vLinks(4*nx*ny,0)
{

/*
 * compute vertex links, and coordinates
 */
	computeVertexLinks();
	computeVertexCoordinates();
	computeLocalVertexCoordinates();
}

void Quadrulator::computeLocalVertexCoordinates(){
	for(size_t e=0;e<numQuad;e++){
		for(size_t i=0;i<npe;i++){
			size_t node = vLinks[e*npe+i];
			localCoord[3*npe*e+3*i]   = coord[3*node];
			localCoord[3*npe*e+3*i+1] = coord[3*node+1];
			localCoord[3*npe*e+3*i+2] = coord[3*node+2];
		}
	}
}

void Quadrulator::computeVertexCoordinates() {
	double x,y;
	size_t p=0;
	for(size_t jy=0;jy<=ny;jy++){
		y = jy * hy + y0;
		for(size_t ix=0;ix<=nx;ix++){
			x = ix * hx + x0;
			coord[p]=x;
			coord[p+1]=y;
			coord[p+2]=0.0;
			p += 3;
		}
	}
}

void Quadrulator::computeVertexLinks(){
	size_t one, two, three, four;
	size_t p = 0;
	for(size_t jy=0;jy<ny;jy++){
		for(size_t ix=0;ix<nx;ix++){
			// quad --> counter clockwise local connectivity
			one = jy * (nx + 1) + ix;
			two = one + 1;
			three = (jy + 1) * (nx+1) + ix + 1;
			four = three - 1;
			vLinks[p]=one;
			vLinks[p+1]=two;
			vLinks[p+2]=three;
			vLinks[p+3]=four;
			p += 4;
		}
	}
}
