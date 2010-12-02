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
  enum FieldType {VOLUME=0, ID, PROC_NUM, WEIGHTED_VOLUME, DILATATION, DAMAGE, E_DP, PLASTIC_CONSISTENCY, NUM_NEIGHBORS, COORDINATES, DISPLACEMENT, CURRENT_COORINATES, VELOCITY, ACCELERATION, FORCE, FORCE_DENSITY, BOND_DAMAGE, DEFAULT_FIELDTYPE};
    enum FieldLength {SCALAR=1, VECTOR2D=2, VECTOR3D=3, BOND=4};
    enum FieldStateArchitecture {STATELESS=0, STATEFUL=1}; 
    enum FieldStep {STEP_N=0, STEP_NP1=1, STEP_NONE=2};

private:
	FieldType type;
	FieldLength length;
    FieldStateArchitecture stateArch;
	std::string label;

public:
	FieldSpec(const FieldSpec& c);
    explicit FieldSpec(FieldType t, FieldLength len, const string& label);
    explicit FieldSpec(FieldType t, FieldLength len, FieldStateArchitecture arch, const string& label);
	bool operator == (const FieldSpec& right) const;
	bool operator != (const FieldSpec& right) const;
	bool operator < (const FieldSpec& right) const;
	FieldType getType() const;
	FieldLength getLength() const;
    FieldStateArchitecture getStateArchitecture() const;

	string getLabel() const;
};

// Scalar FieldSpecs
const FieldSpec DEFAULT_FIELDTYPE(FieldSpec::DEFAULT_FIELDTYPE, FieldSpec::SCALAR, FieldSpec::STATELESS, "Default_FieldType");
const FieldSpec VOLUME(FieldSpec::VOLUME,                       FieldSpec::SCALAR, FieldSpec::STATELESS, "Volume");
const FieldSpec ID(Field_NS::FieldSpec::ID,                     FieldSpec::SCALAR, FieldSpec::STATELESS, "Id");
const FieldSpec PROC_NUM(Field_NS::FieldSpec::PROC_NUM,         FieldSpec::SCALAR, FieldSpec::STATELESS, "Proc_Num");
const FieldSpec DAMAGE(FieldSpec::DAMAGE,                       FieldSpec::SCALAR, FieldSpec::STATEFUL,  "Damage");
const FieldSpec WEIGHTED_VOLUME(FieldSpec::WEIGHTED_VOLUME,     FieldSpec::SCALAR, FieldSpec::STATELESS, "Weighted_Volume");
const FieldSpec DILATATION(FieldSpec::DILATATION,               FieldSpec::SCALAR, FieldSpec::STATEFUL,  "Dilatation");
const FieldSpec NUM_NEIGHBORS(FieldSpec::NUM_NEIGHBORS,         FieldSpec::SCALAR, FieldSpec::STATELESS, "Num_Neighbors");
const FieldSpec LAMBDA(FieldSpec::PLASTIC_CONSISTENCY,          FieldSpec::SCALAR, FieldSpec::STATEFUL,  "Lambda");
const FieldSpec DEVIATORIC_PLASTIC_EXTENSION(FieldSpec::E_DP,   FieldSpec::BOND,   FieldSpec::STATEFUL,  "Deviatoric_Plastic_Extension");

// Vector FieldSpecs
const FieldSpec COORD3D(FieldSpec::COORDINATES,            FieldSpec::VECTOR3D, FieldSpec::STATELESS, "Coordinates");
const FieldSpec DISPL3D(FieldSpec::DISPLACEMENT,           FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Displacement");
const FieldSpec CURCOORD3D(FieldSpec::CURRENT_COORDINATES, FiledSec::VECTOR3D,  FieldSpec::STATEFUL,  "Current_Coordinates");
const FieldSpec ACCEL3D(FieldSpec::ACCELERATION,           FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Acceleration");
const FieldSpec FORCE3D(FieldSpec::FORCE,                  FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Force");
const FieldSpec FORCE_DENSITY3D(FieldSpec::FORCE_DENSITY,  FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Force Density");

// Bond FieldSpecs
const FieldSpec BOND_DAMAGE(FieldSpec::BOND_DAMAGE, FieldSpec::BOND, FieldSpec::STATEFUL, "Bond_Damage");

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
