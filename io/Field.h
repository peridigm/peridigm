/*! \file Field.h */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#ifndef FIELD_H_
#define FIELD_H_
#include "Pd_shared_ptr_Array.h"
#include <stdexcept>
#include <string>
#include <map>
#include <iostream>
#include <sstream>

namespace Field_NS {

using std::tr1::shared_ptr;
using std::string;
using std::map;

/**
 *
 * FieldSpec Declaration
 */

struct FieldSpec {
public:
  enum FieldType {
	  VOLUME=0,
	  ID,
	  PROC_NUM,
	  WEIGHTED_VOLUME,
	  DILATATION,
	  DAMAGE,
	  E_DP,
	  E_DB,
	  PLASTIC_CONSISTENCY,
	  SHEAR_CORRECTION_FACTOR,
	  NUM_NEIGHBORS,
	  COORDINATES,
	  DISPLACEMENT,
	  CURRENT_COORDINATES,
	  VELOCITY,
	  ACCELERATION,
	  FORCE,
	  FORCE_DENSITY,
	  CONTACT_FORCE,
	  CONTACT_FORCE_DENSITY,
	  BOND_DAMAGE,
	  DEFAULT_FIELDTYPE
  };
  enum FieldLength {SCALAR=1, VECTOR2D=2, VECTOR3D=3, BOND=4};
  enum FieldStateArchitecture {STATELESS=0, STATEFUL=1}; 
  enum FieldStep {STEP_N=0, STEP_NP1=1, STEP_NONE=2};

private:
	FieldType type;
	FieldLength length;
        FieldStateArchitecture stateArch;
	std::string label;

public:
	FieldSpec();
	FieldSpec(const FieldSpec& c);
        explicit FieldSpec(FieldType t, FieldLength len, const string& label);
        explicit FieldSpec(FieldType t, FieldLength len, FieldStateArchitecture arch, const string& label);
	bool operator == (const FieldSpec& right) const;
	bool operator != (const FieldSpec& right) const;
	bool operator < (const FieldSpec& right) const;
        std::ostream& print(std::ostream& os) const;
	FieldType getType() const;
	FieldLength getLength() const;
        FieldStateArchitecture getStateArchitecture() const;
	string getLabel() const;
};

inline std::ostream& operator<<(std::ostream& os, const FieldSpec& fs) {
  return fs.print(os);
}

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
const FieldSpec SHEAR_CORRECTION_FACTOR(FieldSpec::SHEAR_CORRECTION_FACTOR, FieldSpec::SCALAR, FieldSpec::STATELESS, "Shear_Correction_Factor");

// Vector FieldSpecs
const FieldSpec COORD3D(FieldSpec::COORDINATES,            FieldSpec::VECTOR3D, FieldSpec::STATELESS, "Coordinates");
const FieldSpec DISPL3D(FieldSpec::DISPLACEMENT,           FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Displacement");
const FieldSpec CURCOORD3D(FieldSpec::CURRENT_COORDINATES, FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Current_Coordinates");
const FieldSpec VELOC3D(FieldSpec::VELOCITY,               FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Velocity");
const FieldSpec ACCEL3D(FieldSpec::ACCELERATION,           FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Acceleration");
const FieldSpec FORCE3D(FieldSpec::FORCE,                  FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Force");
const FieldSpec FORCE_DENSITY3D(FieldSpec::FORCE_DENSITY,  FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Force Density");
const FieldSpec CONTACT_FORCE3D(FieldSpec::CONTACT_FORCE,  FieldSpec::VECTOR3D, FieldSpec::STATEFUL,  "Contact Force");
const FieldSpec CONTACT_FORCE_DENSITY3D(FieldSpec::CONTACT_FORCE_DENSITY, FieldSpec::VECTOR3D, FieldSpec::STATEFUL, "Contact Force Density");

// Bond FieldSpecs
const FieldSpec BOND_DAMAGE(FieldSpec::BOND_DAMAGE,           FieldSpec::BOND, FieldSpec::STATEFUL, "Bond_Damage");
const FieldSpec DEVIATORIC_PLASTIC_EXTENSION(FieldSpec::E_DP, FieldSpec::BOND, FieldSpec::STATEFUL, "Deviatoric_Plastic_Extension");
const FieldSpec DEVIATORIC_BACK_EXTENSION(FieldSpec::E_DB, FieldSpec::BOND, FieldSpec::STATEFUL, "Deviatoric_Back_Extension");

// Map from string ID to FieldSpec, for all fieldspecs
struct FieldSpecMap {
  static std::map<string, FieldSpec> create_map() {
    std::map<string,FieldSpec> mymap;
    mymap[DEFAULT_FIELDTYPE.getLabel()]            = DEFAULT_FIELDTYPE;
    mymap[VOLUME.getLabel()]                       = VOLUME;
    mymap[ID.getLabel()]                           = ID;
    mymap[PROC_NUM.getLabel()]                     = PROC_NUM;
    mymap[DAMAGE.getLabel()]                       = DAMAGE;
    mymap[WEIGHTED_VOLUME.getLabel()]              = WEIGHTED_VOLUME;
    mymap[DILATATION.getLabel()]                   = DILATATION;
    mymap[NUM_NEIGHBORS.getLabel()]                = NUM_NEIGHBORS;
    mymap[LAMBDA.getLabel()]                       = LAMBDA;
    mymap[SHEAR_CORRECTION_FACTOR.getLabel()]      = SHEAR_CORRECTION_FACTOR;
    mymap[COORD3D.getLabel()]                      = COORD3D;
    mymap[DISPL3D.getLabel()]                      = DISPL3D;
    mymap[CURCOORD3D.getLabel()]                   = CURCOORD3D;
    mymap[VELOC3D.getLabel()]                      = VELOC3D;
    mymap[ACCEL3D.getLabel()]                      = ACCEL3D;
    mymap[FORCE3D.getLabel()]                      = FORCE3D;
    mymap[FORCE_DENSITY3D.getLabel()]              = FORCE_DENSITY3D;
    mymap[CONTACT_FORCE3D.getLabel()]              = CONTACT_FORCE3D;
    mymap[CONTACT_FORCE_DENSITY3D.getLabel()]      = CONTACT_FORCE_DENSITY3D;
    mymap[BOND_DAMAGE.getLabel()]                  = BOND_DAMAGE;
    mymap[DEVIATORIC_PLASTIC_EXTENSION.getLabel()] = DEVIATORIC_PLASTIC_EXTENSION;
    mymap[DEVIATORIC_BACK_EXTENSION.getLabel()]    = DEVIATORIC_BACK_EXTENSION ;
    return mymap;
  };
  static const std::map<string, FieldSpec> Map;  
};  

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
