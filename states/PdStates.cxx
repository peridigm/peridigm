/*
 * PdStates.cxx
 *
 *  Created on: Oct 20, 2009
 *      Author: jamitch
 */

#include "CellArray.h"
#include "State.h"
#include "FieldData.h"
#include "FieldSpec.h"
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <cstdio>

using namespace PdStates;
using std::vector;
using std::map;
using std::set;



/**
 * State Implementation; Derived classes must implement operator
 */

/**
 * CellArray Implementation
 */

/**
 * Constructor
 */

CellArray::CellArray(const std::set<FieldSpec>& fieldSpecs, shared_ptr<int>& ids) : vector<int>::iterator(cellIds.get()), cellIds(ids), fieldMap(), numCells(0) {
	// compute numcells
	int* ptr = ids.get();
	numCells = sizeof ptr / sizeof(int);
	// Initialize fieldMap
	for(set<FieldSpec>::iterator iter=fieldSpecs.begin();iter!=fieldSpecs.end();iter++){
		FieldData fd(cellIds,*iter);
		fieldMap.insert(std::pair<FieldSpec,FieldData>(*iter,fd));
	}
}

int CellArray::getNumCells() const { return getNumCells(); }

shared_ptr<int> CellArray::getNeighborhood(int cellId){
	std::cout << "CellArray::getNeighborhood needs to be implemented " << std::endl;
	std::exit(1);
	return shared_ptr<int>();
}

/**
 * FieldData Implementation
 */
FieldData::FieldData(shared_ptr<int>& ids, const FieldSpec& spec) : cellIds(ids), fieldSpec(spec){}


