/*
 * Q2Cylinder.cpp
 *
 *  Created on: Jan 4, 2011
 *      Author: jamitch
 */

Q2Cylinder::Q2Cylinder(const Teuchos::RCP<Teuchos::ParameterList>& peridigmParams)
{
	Teuchos::RCP<Teuchos::ParameterList> discParams = Teuchos::rcp(&(peridigmParams->sublist("Problem").sublist("Discretization")), false);
	  if (discParams->isSublist("TensorProductCylinderMeshGenerator")){
	    Teuchos::RCP<Teuchos::ParameterList> pdQuickGridParamList = Teuchos::rcp(&(params->sublist("TensorProductCylinderMeshGenerator")), false);
	    double innerRadius    = pdQuickGridParamList->get<double>("Inner Radius");
	    double outerRadius    = pdQuickGridParamList->get<double>("Outer Radius");
	    double cylinderLength = pdQuickGridParamList->get<double>("Cylinder Length");
	    int numRings          = pdQuickGridParamList->get<int>("Number Points Radius");
	    double xC             = pdQuickGridParamList->get<double>("Ring Center x");
	    double yC             = pdQuickGridParamList->get<double>("Ring Center y");
	    double zStart         = pdQuickGridParamList->get<double>("Z Origin");

	    // Create 2d Ring
	    std::valarray<double> center(0.0,3);
	    center[0] = xC;
	    center[1] = yC;
	    center[2] = 0;

	  }
	  else { // ERROR
	    TEST_FOR_EXCEPT_MSG(true, "PeridigmNS::InitialConditionsNS::Q2Cylinder-->Invalid Type in PdQuickGridDiscretization ");
	  }


}
