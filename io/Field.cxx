/*! \file Field.cxx */

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
