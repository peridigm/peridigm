/*
 * Peridigm_DecompositionStates.cxx
 *
 *  Created on: Apr 27, 2010
 *      Author: jamitch
 */

#include "Peridigm_DecompositionStates.hpp"
#include "Epetra_MultiVector.h"


PeridigmNS::DecompositionStates::DecompositionStates()
:
scalarStateVarMap(), vectorStateVarMap(), scalarStateBondVarMap(), nextScalarVarIndex(0), nextVectorVarIndex(0), nextScalarBondVarIndex(0)
{

	/**
	 * Scalar state variables
	 */
	scalarStateVarMap["Weighted Volume"]     = nextScalarVarIndex++;
	scalarStateVarMap["Dilatation"]          = nextScalarVarIndex++;
	scalarStateVarMap["Damage"]              = nextScalarVarIndex++;

	/**
	 * Vector state variables
	 */
	vectorStateVarMap["Current Position"]     = nextVectorVarIndex++;


}

const std::string& PeridigmNS::DecompositionStates::getScalarStateName(int pos) const {

	std::map<std::string,int>::const_iterator iter = scalarStateVarMap.begin();
	std::map<std::string,int>::const_iterator end = scalarStateVarMap.end();
	for(;iter!=end;iter++){
		int index = iter->second;
		if(index==pos){
			return iter->first;
		}
	}

	throw std::range_error("PeridigmNS::PNSeridigm_DecompositionStates::getScalarStateName(invalid position).");
}

const std::string& PeridigmNS::DecompositionStates::getVectorStateName(int pos) const {

	std::map<std::string,int>::const_iterator iter = vectorStateVarMap.begin();
	std::map<std::string,int>::const_iterator end = vectorStateVarMap.end();
	for(;iter!=end;iter++){
		int index = iter->second;
		if(index==pos){
			return iter->first;;
		}
	}

	throw std::range_error("PeridigmNS::Peridigm_DecompositionStates::getVectorStateName(invalid position).");
}

const std::string& PeridigmNS::DecompositionStates::getScalarStateBondVarName(int pos) const {

	std::map<std::string,int>::const_iterator iter = scalarStateBondVarMap.begin();
	std::map<std::string,int>::const_iterator end = scalarStateBondVarMap.end();
	for(;iter!=end;iter++){
		int index = iter->second;
		if(index==pos){
			return iter->first;;
		}
	}

	throw std::range_error("PeridigmNS::Peridigm_DecompositionStates::getScalarStateBondVarName(invalid position).");
}

void PeridigmNS::DecompositionStates::addScalarStateVariable(const std::string& name) {
	scalarStateVarMap[name] = nextScalarVarIndex++;
}

void PeridigmNS::DecompositionStates::addVectorStateVariable(const std::string& name) {
	vectorStateVarMap[name] = nextVectorVarIndex++;
}

void PeridigmNS::DecompositionStates::addScalarStateBondVariable(const std::string& name) {
	scalarStateBondVarMap[name] = nextScalarBondVarIndex++;
}

int PeridigmNS::DecompositionStates::getNumScalarStateVariables() const {
	return scalarStateVarMap.size();
}

int PeridigmNS::DecompositionStates::getNumVectorStateVariables() const {
	return vectorStateVarMap.size();
}

int PeridigmNS::DecompositionStates::getNumScalarStateBondVariables() const {
	return scalarStateBondVarMap.size();
}

std::pair<int,double*> PeridigmNS::DecompositionStates::extractStrideView(Epetra_MultiVector& data) const {
	  double* dataView;
	  int dataStride;
	  data.ExtractView(&dataView, &dataStride);
	  return std::make_pair<int,double*>(dataStride,dataView);
}


double* PeridigmNS::DecompositionStates::
extractWeightedVolumeView(std::pair<int,double*>& strideView) const {
	int scalarConstitutiveDataStride=strideView.first;
	double* scalarConstitutiveDataView = strideView.second;
	std::map<std::string,int>::const_iterator iter = scalarStateVarMap.find("Weighted Volume");
	if(iter==scalarStateVarMap.end()){
		throw std::range_error("PeridigmNS::Peridigm_DecompositionStates::extractWeightedVolumeView(...,\"Weighted Volume\") \"Does not exist\"");
		return NULL;
	}
	int index = iter->second;
	double* weightedVolume = scalarConstitutiveDataView + scalarConstitutiveDataStride*index;
	return weightedVolume;
}

double* PeridigmNS::DecompositionStates::
extractDilatationView(std::pair<int,double*>& strideView) const {
	int scalarConstitutiveDataStride=strideView.first;
	double* scalarConstitutiveDataView = strideView.second;
	std::map<std::string,int>::const_iterator iter = scalarStateVarMap.find("Dilatation");
	if(iter==scalarStateVarMap.end()){
		throw std::range_error("PeridigmNS::Peridigm_DecompositionStates::extractDilatationView(...,\"Dilatation\") \"Does not exist\"");
		return NULL;
	}
	int index = iter->second;
	double* dilatation = scalarConstitutiveDataView + scalarConstitutiveDataStride*index;
	return dilatation;
}

double* PeridigmNS::DecompositionStates::
extractDamageView(std::pair<int,double*>& strideView) const {
	int scalarConstitutiveDataStride=strideView.first;
	double* scalarConstitutiveDataView = strideView.second;
	std::map<std::string,int>::const_iterator iter = scalarStateVarMap.find("Damage");
	if(iter==scalarStateVarMap.end()){
		throw std::range_error("PeridigmNS::Peridigm_DecompositionStates::extractDamageView(...,\"Damage\") \"Does not exist\"");
		return NULL;
	}
	int index = iter->second;
	double* dilatation = scalarConstitutiveDataView + scalarConstitutiveDataStride*index;
	return dilatation;
}

double* PeridigmNS::DecompositionStates::
extractCurrentPositionView(std::pair<int,double*>& strideView) const {
	int vectorConstitutiveDataStride = strideView.first;
	double* vectorConstitutiveDataView = strideView.second;
	std::map<std::string,int>::const_iterator iter = vectorStateVarMap.find("Current Position");
	if(iter==vectorStateVarMap.end()){
		throw std::range_error("PeridigmNS::Peridigm_DecompositionStates::extractCurrentPositionView(...,\"Current Position\") \"Does not exist\"");
		return NULL;
	}
	int index = iter->second;
	double* y = vectorConstitutiveDataView + vectorConstitutiveDataStride*index;
	return y;

}

double* PeridigmNS::DecompositionStates::
extractScalarBondVariable(std::pair<int,double*>& strideView, const std::string& name) const {
	int stride = strideView.first;
	double* view = strideView.second;
	std::map<std::string,int>::const_iterator iter = scalarStateBondVarMap.find(name);
	if(iter==scalarStateBondVarMap.end()){
		throw std::range_error("PeridigmNS::Peridigm_DecompositionStates::extractScalarBondVariable"+name+"(...,\"Current Position\") \"Does not exist\"");
		return NULL;
	}
	int index = iter->second;
	double* y = view + stride*index;
	return y;

}
