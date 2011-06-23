/*! \file InitialCondition.cpp */

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

#include "InitialCondition.hpp"
#include "Q2Cylinder.hpp"
#include <map>
#include <string>

namespace InitialConditionsNS {

using std::map;
using std::string;

/*
 * This is a function pointer to a generic builder
 */
typedef RCP<InitialCondition> (*BuilderFunctionPointer)(const Teuchos::ParameterList&);



/*
 * Forward declarations of function pointers to builders
 */
RCP<InitialCondition> getQ2CylinderIC(const Teuchos::ParameterList& peridigmParams);

/*
 * Creates map of support functions
 */
map<string,BuilderFunctionPointer> getSupport() {
	map<string,BuilderFunctionPointer> icSupport;
	icSupport["Q2Cylinder"] = getQ2CylinderIC;
	return icSupport;
}

RCP<InitialCondition> getInstance(const Teuchos::ParameterList& peridigmParams) {
	map<string,BuilderFunctionPointer> icSupport = getSupport();
	const Teuchos::ParameterList& params = peridigmParams.sublist("Problem").sublist("Initial Conditions");
	string icType = params.get<string>("Type");
	BuilderFunctionPointer builder = icSupport[icType];
	return builder(peridigmParams);
}

RCP<InitialCondition> getQ2CylinderIC(const Teuchos::ParameterList& peridigmParams) {
	double vr0, vr1, vz0, z0, a;
	UTILITIES::Vector3D center;

	const Teuchos::ParameterList& params = peridigmParams.sublist("Problem").sublist("Initial Conditions");

	TEST_FOR_EXCEPT_MSG(params.get<string>("Type") != "Q2Cylinder", "PeridigmNS::InitialConditionsNS::getQ2CylinderIC -- invalid \'Type\'");
	if (params.isSublist("IC Params")){
		const Teuchos::ParameterList& ICParams = params.sublist("IC Params");
		TEST_FOR_EXCEPT_MSG(ICParams.get<string>("Geometry") != "TensorProductCylinderMeshGenerator",
				"PeridigmNS::InitialConditionsNS::Q2Cylinder -- Invalid \'Geometry\'");
		vr0 = ICParams.get<double>("v_r0");
		vr1 = ICParams.get<double>("v_r1");
		vz0 = ICParams.get<double>("v_z0");
	}
	else { // ERROR
		TEST_FOR_EXCEPT_MSG(true, "PeridigmNS::InitialConditionsNS::getQ2CylinderIC-->missing \'IC Params\' sublist.");
	}

	const Teuchos::ParameterList& discParams = peridigmParams.sublist("Problem").sublist("Discretization");
	if (discParams.isSublist("TensorProductCylinderMeshGenerator")){
		const Teuchos::ParameterList& pdQuickGridParamList = discParams.sublist("TensorProductCylinderMeshGenerator");
		double cylinderLength = pdQuickGridParamList.get<double>("Cylinder Length");
		double xC             = pdQuickGridParamList.get<double>("Ring Center x");
		double yC             = pdQuickGridParamList.get<double>("Ring Center y");
		double zStart         = pdQuickGridParamList.get<double>("Z Origin");
		/*
		 * Initialize parameters for computing initial conditions
		 */
		z0=zStart;
		a=cylinderLength/2.0;
		center[0]=xC;
		center[1]=yC;
	}
	else { // ERROR
		TEST_FOR_EXCEPT_MSG(true, "PeridigmNS::InitialConditionsNS::Q2Cylinder -- invalid \'Type\' for");
	}
	return RCP<InitialCondition>(new Q2Cylinder(vr0, vr1, vz0, z0, a, center));
}

}  // namespace InitialConditionsNS

