/*
 * Peridigm_DecompositionStates.hpp
 *
 *  Created on: Apr 27, 2010
 *      Author: jamitch
 */


#ifndef PERIDIGM_DECOMPOSITIONSTATES_HPP_
#define PERIDIGM_DECOMPOSITIONSTATES_HPP_
#include <map>
#include <string>
#include <stdexcept>


class Epetra_MultiVector;

namespace Peridigm {

class DecompositionStates {

public:
	DecompositionStates();
	const std::string& getScalarStateName(int pos) const;
	const std::string& getVectorStateName(int pos) const;
	const std::string& getScalarStateBondVarName(int pos) const;
	std::pair<int,double*> extractStrideView(Epetra_MultiVector& data) const;
	double* extractWeightedVolumeView(std::pair<int,double*>& scalarStrideView) const;
	double* extractDilatationView(std::pair<int,double*>& scalarStrideView) const;
	double* extractDamageView(std::pair<int,double*>& scalarStrideView) const;
	double* extractCurrentPositionView(std::pair<int,double*>& vectorStrideView) const;
	double* extractScalarBondVariable(std::pair<int,double*>& vectorStrideView, const std::string& name) const;


	void addScalarStateVariable(const std::string& name);
	void addVectorStateVariable(const std::string& name);
	void addScalarStateBondVariable(const std::string& name);

	int getNumScalarStateVariables() const;
	int getNumVectorStateVariables() const;
	int getNumScalarStateBondVariables() const;

private:
	std::map<std::string,int> scalarStateVarMap, vectorStateVarMap, scalarStateBondVarMap;
	int nextScalarVarIndex, nextVectorVarIndex, nextScalarBondVarIndex;
};

}
#endif /* PERIDIGM_DECOMPOSITIONSTATES_HPP_ */
