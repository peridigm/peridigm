/*
 * FieldData.h
 *
 *  Created on: Oct 20, 2009
 *      Author: jamitch
 */

#ifndef FIELDDATA_H_
#define FIELDDATA_H_

namespace PdStates {

using std::tr1::shared_ptr;
class FieldSpec;

/*
 * Default constructor and destructor
 */
class FieldData {
private:
	shared_ptr<int> cellIds;
	const FieldSpec& fieldSpec;

public:
	FieldData(shared_ptr<int>& cellIds, const FieldSpec& f);

};

}

#endif /* FIELDDATA_H_ */
