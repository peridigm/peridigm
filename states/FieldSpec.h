/*
 * PdField.h
 *
 *  Created on: Oct 19, 2009
 *      Author: jamitch
 */

#ifndef PDFIELD_H_
#define PDFIELD_H_

namespace PdStates {

/**
 *
 * Default constructor, destructor
 */

struct FieldSpec {
public:
	enum FieldType {VOLUME=0, WEIGHTED_VOLUME, COORDINATES, DISPLACEMENT, VELOCITY, ACCELERATION, FORCE};
	enum FieldLength {SCALAR=1, VECTOR2D=2, VECTOR3D=3};

private:
	FieldType type;
	FieldLength length;

public:
	explicit FieldSpec(FieldType t, FieldLength len);
	bool operator == (const FieldSpec& right) const;
	bool operator != (const FieldSpec& right) const;
	bool operator < (const FieldSpec& right) const;
	FieldType getType() const;
	FieldLength getLength() const;
};

const FieldSpec VOLUME(FieldSpec::VOLUME,FieldSpec::SCALAR);
const FieldSpec WEIGHTED_VOLUME(FieldSpec::WEIGHTED_VOLUME,FieldSpec::SCALAR);

const FieldSpec COORD1D(FieldSpec::COORDINATES,FieldSpec::SCALAR);
const FieldSpec COORD2D(FieldSpec::COORDINATES,FieldSpec::VECTOR2D);
const FieldSpec COORD3D(FieldSpec::COORDINATES,FieldSpec::VECTOR3D);

const FieldSpec DISPL1D(FieldSpec::DISPLACEMENT,FieldSpec::SCALAR);
const FieldSpec DISPL2D(FieldSpec::DISPLACEMENT,FieldSpec::VECTOR2D);
const FieldSpec DISPL3D(FieldSpec::DISPLACEMENT,FieldSpec::VECTOR3D);

const FieldSpec VELOC1D(FieldSpec::VELOCITY,FieldSpec::SCALAR);
const FieldSpec VELOC2D(FieldSpec::VELOCITY,FieldSpec::VECTOR2D);
const FieldSpec VELOC3D(FieldSpec::VELOCITY,FieldSpec::VECTOR3D);

const FieldSpec ACCEL1D(FieldSpec::ACCELERATION,FieldSpec::SCALAR);
const FieldSpec ACCEL2D(FieldSpec::ACCELERATION,FieldSpec::VECTOR2D);
const FieldSpec ACCEL3D(FieldSpec::ACCELERATION,FieldSpec::VECTOR3D);

const FieldSpec FORCE1D(FieldSpec::FORCE,FieldSpec::SCALAR);
const FieldSpec FORCE2D(FieldSpec::FORCE,FieldSpec::VECTOR2D);
const FieldSpec FORCE3D(FieldSpec::FORCE,FieldSpec::VECTOR3D);

}

#endif /* PDFIELD_H_ */
