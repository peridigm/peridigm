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
#include <stdexcept>
#include <string>
#include <map>
#include "utilities/Array.h"

using std::tr1::shared_ptr;
using UTILITIES::Array;
using std::string;
using std::size_t;
using std::map;

namespace Field_ENUM {

enum Relation {
  NODE=0,
  ELEMENT,
  BOND,
  BLOCK,
  GLOBAL,
  RELATION_UNDEFINED
};

enum Type {
  VOLUME=0,
  DENSITY,
  GID,
  BLOCK_ID,
  PROC_NUM,
  WEIGHTED_VOLUME,
  RADIUS,
  NEIGHBORHOOD_VOLUME,
  NUMBER_OF_NEIGHBORS,
  CRITICAL_TIME_STEP,
  DILATATION,
  DAMAGE,
  E_DP,
  E_DB,
  PLASTIC_CONSISTENCY,
  NORM_DEVIATORIC_FORCE_STATE,
  SHEAR_CORRECTION_FACTOR,
  NUM_NEIGHBORS,
  COORDINATES,
  DISPLACEMENT,
  CURRENT_COORDINATES,
  VELOCITY,
  ACCELERATION,
  BC_MASK,
  FORCE,
  FORCE_DENSITY,
  CONTACT_FORCE,
  CONTACT_FORCE_DENSITY,
  RESIDUAL,
  BOND_DAMAGE,
  PARTIAL_VOLUME,
  TYPE_UNDEFINED,
  ANGULAR_MOMENTUM,
  LINEAR_MOMENTUM,
  KINETIC_ENERGY,
  STRAIN_ENERGY,
  STRAIN_ENERGY_DENSITY	
};

enum Length {
  LENGTH_UNDEFINED=0,
  SCALAR=1,
  VECTOR2D=2,
  VECTOR3D=3
};

enum Temporal {
  SCRATCH=0,
  CONSTANT,
  TWO_STEP,
  TEMPORAL_UNDEFINED
};

enum ParallelTopology {
  OWNED=0,
  OVERLAP,
  NONE,
  PARALLEL_TOPOLOGY_UNDEFINED
};

enum Step {
  STEP_N=0,
  STEP_NP1=1,
  STEP_NONE=2,
  STEP_UNDEFINED
};

}

namespace Field_NS {

struct FieldSpec  {
private:
  Field_ENUM::Type type;
  Field_ENUM::Relation relation;
  Field_ENUM::Length length;
  Field_ENUM::ParallelTopology par_top;
  Field_ENUM::Temporal temporal;
  unsigned int id;
  string label;
  FieldSpec(Field_ENUM::Type t, Field_ENUM::Relation r,  Field_ENUM::Length len, Field_ENUM::ParallelTopology p, Field_ENUM::Temporal temp, const string& label);

public:
  FieldSpec();
  FieldSpec(const FieldSpec& c);
  explicit FieldSpec(Field_ENUM::Type t, Field_ENUM::Relation r,  Field_ENUM::Length len, Field_ENUM::Temporal temp, const string& label);
  const FieldSpec& operator=(const FieldSpec& c);
  bool operator == (const FieldSpec& right) const;
  bool operator != (const FieldSpec& right) const;
  bool operator < (const FieldSpec& right) const;
  std::ostream& print(std::ostream& os) const;
  const Field_ENUM::Type getType() const { return type; }
  const Field_ENUM::Relation getRelation() const { return relation; }
  const Field_ENUM::Length getLength() const { return length; }
  const Field_ENUM::ParallelTopology get_parallel_topology() const { return par_top; }
  const Field_ENUM::Temporal get_temporal() const { return temporal; }
  const unsigned int getId() const { return id; }
  const string getLabel() const { return label; }
  const FieldSpec get_overlap_spec() const;
  const FieldSpec get_override(Field_ENUM::Temporal temp) const;
};

inline std::ostream& operator<<(std::ostream& os, const FieldSpec& fs) {
  return fs.print(os);
}

/*
 * UNDEFINED FIELDSPEC
 */
const FieldSpec FIELDSPEC_UNDEFINED(Field_ENUM::TYPE_UNDEFINED,Field_ENUM::RELATION_UNDEFINED,Field_ENUM::LENGTH_UNDEFINED,Field_ENUM::TEMPORAL_UNDEFINED,"FIELDSPEC_UNDEFINED");

template<typename T>
class Field : public FieldSpec, public Array<T> {
private:
  size_t num_points;

public:
  Field() : FieldSpec(FIELDSPEC_UNDEFINED), Array<T>(0*FIELDSPEC_UNDEFINED.getLength()), num_points(0) {}
  Field(const FieldSpec& c, size_t numPoints) : FieldSpec(c), Array<T>(numPoints*c.getLength()), num_points(numPoints) {}
  Field(const FieldSpec& c, shared_ptr<T> data, size_t numPoints) :
  FieldSpec(c), Array<T>(numPoints*c.getLength(),data), num_points(numPoints) {}
  size_t get_num_points() const { return num_points; }
};


template<typename T>
class TemporalField {

public:
  TemporalField(const FieldSpec& c=FIELDSPEC_UNDEFINED,  size_t numPoints=0): num_points(numPoints), N(c,numPoints), NP1(c,numPoints) {}
  size_t get_num_points() const { return num_points; }
  Field<T> getField(Field_ENUM::Step step) const {
    switch(step){
      case Field_ENUM::STEP_N:
        return N;
        break;
      case Field_ENUM::STEP_NP1:
        return NP1;
        break;
      default:
        std::string m = "TemporalField::TemporalField(Invalid value for argument Step)\n";
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
  size_t num_points;
  Field<T> N, NP1;
};

/*
 * Constant FieldSpecs
 */

/*
 * GLOBAL SCALAR FieldSpecs (scalar fields defined over entire simulation)
 */
const Field_NS::FieldSpec GLOBAL_KINETIC_ENERGY(Field_ENUM::KINETIC_ENERGY,               Field_ENUM::GLOBAL, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Global_Kinetic_Energy");
const Field_NS::FieldSpec GLOBAL_STRAIN_ENERGY(Field_ENUM::STRAIN_ENERGY,                 Field_ENUM::GLOBAL, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Global_Strain_Energy");
const Field_NS::FieldSpec GLOBAL_STRAIN_ENERGY_DENSITY(Field_ENUM::STRAIN_ENERGY_DENSITY, Field_ENUM::GLOBAL, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Global_Strain_Energy_Density");
const Field_NS::FieldSpec GLOBAL_LINEAR_MOMENTUM(Field_ENUM::LINEAR_MOMENTUM,             Field_ENUM::GLOBAL, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Global_Linear_Momentum");
const Field_NS::FieldSpec GLOBAL_ANGULAR_MOMENTUM(Field_ENUM::ANGULAR_MOMENTUM,           Field_ENUM::GLOBAL, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Global_Angular_Momentum");

/*
 * ELEMENT SCALAR FieldSpecs (scalar fields defined over elements)
 */
const Field_NS::FieldSpec VOLUME                  (Field_ENUM::VOLUME,                       Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Volume");
const Field_NS::FieldSpec DENSITY                 (Field_ENUM::DENSITY,                      Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Density");
const Field_NS::FieldSpec GID                     (Field_ENUM::GID,                          Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "ID");
const Field_NS::FieldSpec BLOCK_ID                (Field_ENUM::BLOCK_ID,                     Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "BLOCK_ID");
const Field_NS::FieldSpec PROC_NUM                (Field_ENUM::PROC_NUM,                     Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Proc_Num");
const Field_NS::FieldSpec WEIGHTED_VOLUME         (Field_ENUM::WEIGHTED_VOLUME,              Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Weighted_Volume");
const Field_NS::FieldSpec RADIUS                  (Field_ENUM::RADIUS,                       Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Radius");
const Field_NS::FieldSpec NEIGHBORHOOD_VOLUME     (Field_ENUM::NEIGHBORHOOD_VOLUME,          Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Neighborhood_Volume");
const Field_NS::FieldSpec NUMBER_OF_NEIGHBORS     (Field_ENUM::NUMBER_OF_NEIGHBORS,          Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Number_Of_Neighbors");
const Field_NS::FieldSpec CRITICAL_TIME_STEP      (Field_ENUM::CRITICAL_TIME_STEP,           Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Critical_Time_Step");
const Field_NS::FieldSpec DILATATION              (Field_ENUM::DILATATION,                   Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::TWO_STEP, "Dilatation");
const Field_NS::FieldSpec NUM_NEIGHBORS           (Field_ENUM::NUM_NEIGHBORS,                Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Num_Neighbors");
const Field_NS::FieldSpec LAMBDA                  (Field_ENUM::PLASTIC_CONSISTENCY,          Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::TWO_STEP, "Lambda");
const Field_NS::FieldSpec NORM_TD                 (Field_ENUM::NORM_DEVIATORIC_FORCE_STATE,  Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::SCRATCH,  "Norm_td");
const Field_NS::FieldSpec SHEAR_CORRECTION_FACTOR (Field_ENUM::SHEAR_CORRECTION_FACTOR,      Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "SHEAR_CORRECTION_FACTOR");
const Field_NS::FieldSpec BC_MASK                 (Field_ENUM::BC_MASK,                      Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "BC_MASK");
const Field_NS::FieldSpec DAMAGE                  (Field_ENUM::DAMAGE,                       Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::TWO_STEP, "Damage");
const Field_NS::FieldSpec KINETIC_ENERGY          (Field_ENUM::KINETIC_ENERGY,               Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::TWO_STEP, "Kinetic_Energy");
const Field_NS::FieldSpec STRAIN_ENERGY           (Field_ENUM::STRAIN_ENERGY,                Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::TWO_STEP, "Strain_Energy");
const Field_NS::FieldSpec STRAIN_ENERGY_DENSITY   (Field_ENUM::STRAIN_ENERGY_DENSITY,        Field_ENUM::ELEMENT, Field_ENUM::SCALAR, Field_ENUM::TWO_STEP, "Strain_Energy_Density");

/*
 * NODAL VECTOR3D FieldSpecs (vector fields defined at nodes)
 */
const Field_NS::FieldSpec COORD3D                (Field_ENUM::COORDINATES,           Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::CONSTANT, "Coordinates");
const Field_NS::FieldSpec DISPL3D                (Field_ENUM::DISPLACEMENT,          Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Displacement");
const Field_NS::FieldSpec CURCOORD3D             (Field_ENUM::CURRENT_COORDINATES,   Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Current_Coordinates");
const Field_NS::FieldSpec VELOC3D                (Field_ENUM::VELOCITY,              Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Velocity");
const Field_NS::FieldSpec ACCEL3D                (Field_ENUM::ACCELERATION,          Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Acceleration");
const Field_NS::FieldSpec FORCE3D                (Field_ENUM::FORCE,                 Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Force");
const Field_NS::FieldSpec FORCE_DENSITY3D        (Field_ENUM::FORCE_DENSITY,         Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Force_Density");
const Field_NS::FieldSpec CONTACT_FORCE3D        (Field_ENUM::CONTACT_FORCE,         Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Contact_Force");
const Field_NS::FieldSpec CONTACT_FORCE_DENSITY3D(Field_ENUM::CONTACT_FORCE_DENSITY, Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Contact_Force_Density");
const Field_NS::FieldSpec RESID3D                (Field_ENUM::RESIDUAL,              Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::SCRATCH,  "Residual");
const Field_NS::FieldSpec ANGULAR_MOMENTUM3D     (Field_ENUM::ANGULAR_MOMENTUM,      Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Angular_Momentum");
const Field_NS::FieldSpec LINEAR_MOMENTUM3D      (Field_ENUM::LINEAR_MOMENTUM,       Field_ENUM::NODE, Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP, "Linear_Momentum");

/*
 * BOND SCALAR FieldSpecs
 */
const Field_NS::FieldSpec BOND_DAMAGE                 (Field_ENUM::BOND_DAMAGE,    Field_ENUM::BOND, Field_ENUM::SCALAR, Field_ENUM::TWO_STEP, "Bond Damage");
const Field_NS::FieldSpec PARTIAL_VOLUME              (Field_ENUM::PARTIAL_VOLUME, Field_ENUM::BOND, Field_ENUM::SCALAR, Field_ENUM::CONSTANT, "Partial Volume");
const Field_NS::FieldSpec DEVIATORIC_PLASTIC_EXTENSION(Field_ENUM::E_DP,           Field_ENUM::BOND, Field_ENUM::SCALAR, Field_ENUM::TWO_STEP, "Deviatoric_Plastic_Extension");
const Field_NS::FieldSpec DEVIATORIC_BACK_EXTENSION   (Field_ENUM::E_DB,           Field_ENUM::BOND, Field_ENUM::SCALAR, Field_ENUM::TWO_STEP, "Deviatoric_Back_Extension");

struct FieldSpecMap {
	static std::map<string, FieldSpec> create_map() {
		std::map<string,FieldSpec> mymap;
                // Undefined fieldspec
		mymap[FIELDSPEC_UNDEFINED.getLabel()]          = FIELDSPEC_UNDEFINED;
                // global scalar fieldspecs
		mymap[GLOBAL_KINETIC_ENERGY.getLabel()]        = GLOBAL_KINETIC_ENERGY;
                mymap[GLOBAL_STRAIN_ENERGY.getLabel()]         = GLOBAL_STRAIN_ENERGY;
		mymap[GLOBAL_STRAIN_ENERGY_DENSITY.getLabel()] = GLOBAL_STRAIN_ENERGY_DENSITY;
		mymap[GLOBAL_LINEAR_MOMENTUM.getLabel()]       = GLOBAL_LINEAR_MOMENTUM;
		mymap[GLOBAL_ANGULAR_MOMENTUM.getLabel()]      = GLOBAL_ANGULAR_MOMENTUM;
		// point scalar fieldspecs
		mymap[VOLUME.getLabel()]                       = VOLUME;
		mymap[DENSITY.getLabel()]                      = DENSITY;
		mymap[GID.getLabel()]                          = GID;
		mymap[BLOCK_ID.getLabel()]                     = BLOCK_ID;
		mymap[PROC_NUM.getLabel()]                     = PROC_NUM;
		mymap[WEIGHTED_VOLUME.getLabel()]              = WEIGHTED_VOLUME;
		mymap[RADIUS.getLabel()]                       = RADIUS;
		mymap[NEIGHBORHOOD_VOLUME.getLabel()]          = NEIGHBORHOOD_VOLUME;
		mymap[NUMBER_OF_NEIGHBORS.getLabel()]          = NUMBER_OF_NEIGHBORS;
		mymap[CRITICAL_TIME_STEP.getLabel()]           = CRITICAL_TIME_STEP;
		mymap[DILATATION.getLabel()]                   = DILATATION;
		mymap[NUM_NEIGHBORS.getLabel()]                = NUM_NEIGHBORS;
		mymap[LAMBDA.getLabel()]                       = LAMBDA;
		mymap[NORM_TD.getLabel()]                      = NORM_TD;
		mymap[SHEAR_CORRECTION_FACTOR.getLabel()]      = SHEAR_CORRECTION_FACTOR;
		mymap[BC_MASK.getLabel()]                      = BC_MASK;
		mymap[DAMAGE.getLabel()]                       = DAMAGE;
		mymap[KINETIC_ENERGY.getLabel()]               = KINETIC_ENERGY;
		mymap[STRAIN_ENERGY.getLabel()]                = STRAIN_ENERGY;
		mymap[STRAIN_ENERGY_DENSITY.getLabel()]        = STRAIN_ENERGY_DENSITY;
                // point vector fieldspecs
		mymap[COORD3D.getLabel()]                      = COORD3D;
		mymap[DISPL3D.getLabel()]                      = DISPL3D;
		mymap[CURCOORD3D.getLabel()]                   = CURCOORD3D;
		mymap[VELOC3D.getLabel()]                      = VELOC3D;
		mymap[ACCEL3D.getLabel()]                      = ACCEL3D;
		mymap[FORCE3D.getLabel()]                      = FORCE3D;
		mymap[FORCE_DENSITY3D.getLabel()]              = FORCE_DENSITY3D;
		mymap[CONTACT_FORCE3D.getLabel()]              = CONTACT_FORCE3D;
		mymap[CONTACT_FORCE_DENSITY3D.getLabel()]      = CONTACT_FORCE_DENSITY3D;
		mymap[RESID3D.getLabel()]                      = RESID3D;
		mymap[ANGULAR_MOMENTUM3D.getLabel()]           = ANGULAR_MOMENTUM3D;
		mymap[LINEAR_MOMENTUM3D.getLabel()]            = LINEAR_MOMENTUM3D;
                // bond scalar fieldspecs
		mymap[BOND_DAMAGE.getLabel()]                  = BOND_DAMAGE;
		mymap[PARTIAL_VOLUME.getLabel()]               = PARTIAL_VOLUME;
		mymap[DEVIATORIC_PLASTIC_EXTENSION.getLabel()] = DEVIATORIC_PLASTIC_EXTENSION;
		mymap[DEVIATORIC_BACK_EXTENSION.getLabel()]    = DEVIATORIC_BACK_EXTENSION ;
		return mymap;
	};
	static const std::map<string, FieldSpec> Map;
};

} // Field_NS

#endif // FIELD_H_
