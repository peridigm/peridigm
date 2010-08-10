/*
 * CellArray.h
 *
 *  Created on: Oct 20, 2009
 *      Author: jamitch
 */

#ifndef CELLARRAY_H_
#define CELLARRAY_H_
#include <map>
#include <vector>
#include <set>
#include <iterator>
#include <tr1/memory>

namespace PdStates {

class FieldData;
class FieldSpec;
using std::vector;
using std::map;
using std::set;
using std::iterator;
using std::tr1::shared_ptr;


/**
 * Holds a specification for all fields stored on cells;
 *
 */


class CellArray : public iterator {
private:
	shared_ptr<int> cellIds;
	map<FieldSpec, FieldData> fieldMap;
	int numCells;
public:
	/**
	 * Constructs master set of cells
	 */
	CellArray(const set<FieldSpec>& fieldSpecs, shared_ptr<int>& cellIds);

	/**
	 * Construct a neighborhood around a "this" cell
	 */
	CellArray getNeighborhood();

	int getNumCells() const;

	/**
	 * Fields
	 */
	FieldData& getField(const FieldSpec& spec);
};

}

#endif /* CELLARRAY_H_ */
