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

#include "Field.h"

namespace Field_NS {


/*
 * Static variable that must be initialized
 */
const map<string, FieldSpec> FieldSpecMap::Map =  FieldSpecMap::create_map();

std::ostream& FieldSpec::print(std::ostream& os) const {
	os << label;
	return os;
}

/**
 * FieldSpec Implementation
 */

// Default constructor
FieldSpec::FieldSpec()
:
		type(FIELDSPEC_UNDEFINED.getType()),
		relation(FIELDSPEC_UNDEFINED.getRelation()),
		length(FIELDSPEC_UNDEFINED.getLength()),
		par_top(FIELDSPEC_UNDEFINED.get_parallel_topology()),
		temporal(FIELDSPEC_UNDEFINED.get_temporal()),
		id(FIELDSPEC_UNDEFINED.getId()),
		label(FIELDSPEC_UNDEFINED.getLabel())
{}

FieldSpec::
FieldSpec(Field_ENUM::Type t, Field_ENUM::Relation r,  Field_ENUM::Length len, Field_ENUM::Temporal temp, const string& label)
:
		type(t),    relation(r),        length(len),             par_top(Field_ENUM::OWNED), temporal(temp),
		id(type | (relation<<8) | (length << (8+4)) | (par_top << (8+4+4)) |            (temporal << (8+4+4+4))),
		label(label)
{}

FieldSpec::
FieldSpec(Field_ENUM::Type t, Field_ENUM::Relation r,  Field_ENUM::Length len, Field_ENUM::ParallelTopology p, Field_ENUM::Temporal temp, const string& label)
:
	type(t),    relation(r),        length(len),             par_top(p), temporal(temp),
	id(type | (relation<<8) | (length << (8+4)) | (par_top << (8+4+4)) |            (temporal << (8+4+4+4))),
	label(label)
{}

FieldSpec::FieldSpec(const FieldSpec& c)
:
		type(c.getType()),
		relation(c.getRelation()),
		length(c.getLength()),
		par_top(c.get_parallel_topology()),
		temporal(c.get_temporal()),
		id(c.getId()),
		label(c.getLabel()) {}

const FieldSpec& FieldSpec::operator=(const FieldSpec& right)  {
	if(this != &right){
		type = right.getType();
		relation=right.getRelation();
		length = right.getLength();
		par_top=right.get_parallel_topology();
		temporal = right.get_temporal();
		id=right.getId();
		label = right.getLabel();
	}
	return *this;
}

bool FieldSpec::operator == (const FieldSpec& right) const {
  return (id == right.getId()) ? true : false;
}

bool FieldSpec::operator != (const FieldSpec& right) const {
	return !(*this==right);
}

bool FieldSpec::operator < (const FieldSpec& right) const {
	return id < right.getId();
}

const FieldSpec FieldSpec::get_overlap_spec() const {
	return FieldSpec(type,relation,length,Field_ENUM::OVERLAP,temporal,label);
}

const FieldSpec FieldSpec::get_override(Field_ENUM::Temporal override) const {
	return FieldSpec(type,relation,length,par_top,override,label);
}

} // Field

