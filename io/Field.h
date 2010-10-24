/*
 * Field.h
 *
 *  Created on: Mar 18, 2010
 *      Author: jamitch
 */

#ifndef FIELD_H_
#define FIELD_H_
#include "Pd_shared_ptr_Array.h"
#include <stdexcept>
#include <string>

namespace Field_NS {

using std::tr1::shared_ptr;
using std::string;
/**
 *
 * FieldSpec Declaration
 */

struct FieldSpec {
public:
	enum FieldType {VOLUME=0, ID, PROC_NUM, WEIGHTED_VOLUME, DILATATION, DAMAGE, E_DP, NUM_NEIGHBORS, COORDINATES, DISPLACEMENT, VELOCITY, ACCELERATION, FORCE, DEFAULT_FIELDTYPE};
	enum FieldLength {SCALAR=1, VECTOR2D=2, VECTOR3D=3};
    enum FieldStep {STEP_N=0, STEP_NP1=1, STEP_NONE=2};

private:
	FieldType type;
	FieldLength length;
	std::string label;

public:
	FieldSpec(const FieldSpec& c);
	explicit FieldSpec(FieldType t, FieldLength len, const string& label);
	bool operator == (const FieldSpec& right) const;
	bool operator != (const FieldSpec& right) const;
	bool operator < (const FieldSpec& right) const;
	FieldType getType() const;
	FieldLength getLength() const;

	string getLabel() const;
};

// Scalar FieldSpecs
const FieldSpec DEFAULT_FIELDTYPE(FieldSpec::DEFAULT_FIELDTYPE,FieldSpec::SCALAR, "Default_FieldType");
const FieldSpec VOLUME(FieldSpec::VOLUME,FieldSpec::SCALAR, "Volume");
const FieldSpec ID(Field_NS::FieldSpec::ID,Field_NS::FieldSpec::SCALAR, "Id");
const FieldSpec PROC_NUM(Field_NS::FieldSpec::PROC_NUM,Field_NS::FieldSpec::SCALAR, "Proc_Num");
const FieldSpec DAMAGE(FieldSpec::VOLUME,FieldSpec::SCALAR, "Damage");
const FieldSpec WEIGHTED_VOLUME(FieldSpec::WEIGHTED_VOLUME,FieldSpec::SCALAR, "Weighted_Volume");
const FieldSpec DILATATION(FieldSpec::DILATATION,FieldSpec::SCALAR, "Dilatation");
const FieldSpec NUM_NEIGHBORS(FieldSpec::NUM_NEIGHBORS,FieldSpec::SCALAR, "Num_Neighbors");
const FieldSpec DEVIATORIC_PLASTIC_EXTENSION(FieldSpec::E_DP,FieldSpec::SCALAR, "Deviatoric_Plastic_Extension");

// Vector FieldSpecs
const FieldSpec COORD3D(FieldSpec::COORDINATES,FieldSpec::VECTOR3D, "Coordinates");
const FieldSpec DISPL3D(FieldSpec::DISPLACEMENT,FieldSpec::VECTOR3D, "Displacement");
const FieldSpec VELOC3D(FieldSpec::VELOCITY,    FieldSpec::VECTOR3D, "Velocity");
const FieldSpec ACCEL3D(FieldSpec::ACCELERATION,FieldSpec::VECTOR3D, "Acceleration");
const FieldSpec FORCE3D(FieldSpec::FORCE,FieldSpec::VECTOR3D, "Force");

template<typename T>
class Field : public FieldSpec {

private:
	Pd_shared_ptr_Array<T> data;
	std::size_t numPoints;

public:
    Field(const FieldSpec& c, std::size_t numPoints) : FieldSpec(c), data(numPoints*c.getLength()), numPoints(numPoints){}
    Field(const FieldSpec& c, shared_ptr<T> data, std::size_t numPoints) : FieldSpec(c), data(numPoints*c.getLength(),data), numPoints(numPoints) {}
    std::size_t getNumPoints() const { return numPoints; }
    Pd_shared_ptr_Array<T>& getArray() { return data; }
    const Pd_shared_ptr_Array<T>& getArray() const { return data; }
    void setValue(T value) {
    	const T* const end = data.end();
    	T* f = data.get();
    	for(;f!=end;f++)
    		*f=value;
    }
};


template<typename T>
class TemporalField {
public:

	TemporalField(const FieldSpec& c=DEFAULT_FIELDTYPE, std::size_t numPoints=0): numPoints(numPoints), N(c,numPoints), NP1(c,numPoints) {}
	std::size_t getNumPoints() const { return numPoints; }
    Field<T> getField(FieldSpec::FieldStep step) const {
    	switch(step){
    	case FieldSpec::STEP_N:
    		return N;
    		break;
    	case FieldSpec::STEP_NP1:
    		return NP1;
    		break;
    	default:
    		std::string m = "TemporalField::TemporalField(Invalid value for argument FieldStep)\n";
    		throw std::invalid_argument(m);
    	}
    }

    /*
     * This function swaps the data between N and NP1;
     * Calculations are generally performed on NP1.
     * During time integration, NP1 can be thought of
     * as a scratch array;  At some point, NP1 will satisfy
     * iteration requirements, then, this call
     * advances the time step in order to save NP1 to N.
     * Then the process begins anew.
     */
    void advanceStep() { std::swap(N,NP1); }

private:
	int numPoints;
	Field<T> N, NP1;
};

Field<double>         getVOLUME(shared_ptr<double> data, int numPoints);
Field<double>         getWEIGHTED_VOLUME(int numPoints);
Field<int>            getNUM_NEIGHBORS(int numPoints);
TemporalField<double> getDISPL3D(int numPoints);
TemporalField<double> getVELOC3D(int numPoints);
Field<double>         getFORCE3D(int numPoints);
TemporalField<double> getDILATATION(int numPoints);

}  // namespace Field_NS

#endif /* FIELD_H_ */
