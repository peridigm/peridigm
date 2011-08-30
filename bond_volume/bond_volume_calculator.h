/*
 * bond_volume_calculator.h
 *
 *  Created on: Jun 21, 2011
 *      Author: jamitch
 */

#ifndef BOND_VOLUME_CALCULATOR_H_
#define BOND_VOLUME_CALCULATOR_H_

namespace MATERIAL_EVALUATION {

#include <functional>
using std::binary_function;

struct Bond_Volume_Calculator : public std::binary_function<const double*, const double*, double> {
	virtual ~Bond_Volume_Calculator(){}
	virtual double operator() (const double* P, const double* Q) const = 0;
	virtual double get_horizon() const = 0;
	virtual double get_cell_diagonal() const = 0;
};
}

#endif /* BOND_VOLUME_CALCULATOR_H_ */
