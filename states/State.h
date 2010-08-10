

#ifndef PDSTATE_H_
#define PDSTATE_H_


namespace PdStates {

class FieldData;
class CellArray;
using std::unary_function;

/**
 * This class represents a state by:
 * ** operating on each cell to produce FieldData  -- this can be scalar or vector
 */
class State : public std::unary_function<CellArray&, FieldData&> {
public:
	virtual ~State();
	virtual FieldData& operator()(CellArray& cellArray) const = 0;
};


}
#endif // PDSTATE_H_
