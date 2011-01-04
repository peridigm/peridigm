/*
 * Q2Cylinder.cpp
 *
 *  Created on: Jan 4, 2011
 *      Author: jamitch
 */

#include "Q2Cylinder.hpp"
#include <string>

using std::string;
namespace PeridigmNS {

namespace InitialConditionsNS {


Q2Cylinder::Q2Cylinder(const Teuchos::RCP<Teuchos::ParameterList>& peridigmParams):
		z0(0.0), a(0.0), vr0(0.0), vr1(0.0), vz0(0.0), center()
{

	Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(&(peridigmParams->sublist("Problem").sublist("Initial Conditions")), false);
	TEST_FOR_EXCEPT_MSG(params->get<string>("Type") != "Q2Cylinder", "PeridigmNS::InitialConditionsNS::Q2Cylinder -- invalid \'Type\'");
	if (params->isSublist("IC Params")){
		Teuchos::RCP<Teuchos::ParameterList> ICParams = Teuchos::rcp(&(params->sublist("IC Params")), false);
		TEST_FOR_EXCEPT_MSG(ICParams->get<string>("Geometry") != "TensorProductCylinderMeshGenerator",
				"PeridigmNS::InitialConditionsNS::Q2Cylinder -- Invalid \'Geometry\'");
		vr0 = ICParams->get<double>("v_r0");
		vr1 = ICParams->get<double>("v_r1");
		vz0 = ICParams->get<double>("v_z0");
	}
	else { // ERROR
		TEST_FOR_EXCEPT_MSG(true, "PeridigmNS::InitialConditionsNS::Q2Cylinder-->missing \'IC Params\' sublist.");
	}
	Teuchos::RCP<Teuchos::ParameterList> discParams = Teuchos::rcp(&(peridigmParams->sublist("Problem").sublist("Discretization")), false);
	if (discParams->isSublist("TensorProductCylinderMeshGenerator")){
		Teuchos::RCP<Teuchos::ParameterList> pdQuickGridParamList = Teuchos::rcp(&(discParams->sublist("TensorProductCylinderMeshGenerator")), false);
		double cylinderLength = pdQuickGridParamList->get<double>("Cylinder Length");
		double xC             = pdQuickGridParamList->get<double>("Ring Center x");
		double yC             = pdQuickGridParamList->get<double>("Ring Center y");
		double zStart         = pdQuickGridParamList->get<double>("Z Origin");

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

}

}


}
