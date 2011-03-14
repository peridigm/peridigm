/*
 * Field.cxx
 *
 *  Created on: Mar 18, 2010
 *      Author: jamitch
 */
#include "Field.h"

namespace Field_NS {

/**
 * FieldSpec Implementation
 */
FieldSpec::FieldSpec(FieldType t, FieldLength len, const string& label) : type(t), length(len), stateArch(STATELESS), label(label){}
FieldSpec::FieldSpec(FieldType t, FieldLength len, FieldStateArchitecture arch, const string& label) : type(t), length(len), stateArch(arch), label(label){}
FieldSpec::FieldSpec(const FieldSpec& c) : type(c.getType()), length(c.getLength()), stateArch(c.getStateArchitecture()), label(c.getLabel()) {}
bool FieldSpec::operator == (const FieldSpec& right) const {
  return (type == right.getType() && length == right.getLength() && stateArch == right.getStateArchitecture()) ? true : false;
}
bool FieldSpec::operator != (const FieldSpec& right) const {
	return !(*this==right);
}

bool FieldSpec::operator < (const FieldSpec& right) const {
	if (length < right.getLength())
		return true;
	else if(length == right.getLength())
		return type < right.getType();
	else
		return false;
}
FieldSpec::FieldType FieldSpec::getType() const { return type; }
FieldSpec::FieldLength FieldSpec::getLength() const { return length; }
FieldSpec::FieldStateArchitecture FieldSpec::getStateArchitecture() const { return stateArch; }
string FieldSpec::getLabel() const { return label; }

Field<double> getVOLUME(shared_ptr<double> data, int numPoints) {
    return Field<double>(VOLUME,data,numPoints);
}

Field<double> getWEIGHTED_VOLUME(int numPoints){
	return Field<double>(WEIGHTED_VOLUME,numPoints);
}

Field<int> getNUM_NEIGHBORS(int numPoints){
	return Field<int>(NUM_NEIGHBORS,numPoints);
}

TemporalField<double> getDISPL3D(int numPoints){
	return TemporalField<double>(DISPL3D,numPoints);
}

TemporalField<double> getVELOC3D(int numPoints) {
	return TemporalField<double>(VELOC3D,numPoints);
}

Field<double> getFORCE3D(int numPoints){
	return Field<double>(FORCE3D,numPoints);
}

TemporalField<double> getDILATATION(int numPoints){
	return TemporalField<double>(DILATATION,numPoints);
}

} // Field_NS
