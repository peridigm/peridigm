/*! \file Peridigm_Compute_Error.cpp */

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
// Michael L. Parks      parks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include <vector>

#include "Peridigm_Compute_Error.hpp"
#include "Peridigm_Field.hpp"
#include <sstream>

//! Standard constructor.
PeridigmNS::Compute_Error::Compute_Error(Teuchos::RCP<const Teuchos::ParameterList> params,
                                         Teuchos::RCP<const Epetra_Comm> epetraComm_,
                                         Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_)
  : Compute(params, epetraComm_, computeClassGlobalData_), m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_errorFieldId(-1), m_globalErrorFieldId(-1),
    m_exodusNode1FieldId(-1), m_exodusNode2FieldId(-1), m_exodusNode3FieldId(-1), m_exodusNode4FieldId(-1),
    m_exodusNode5FieldId(-1), m_exodusNode6FieldId(-1), m_exodusNode7FieldId(-1), m_exodusNode8FieldId(-1)
{
  FieldManager& fieldManager = FieldManager::self();
  m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  m_coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  m_errorFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Error");
  m_globalErrorFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Global_Error");
  m_exodusNode1FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Exodus_Node_1");
  m_exodusNode2FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Exodus_Node_2");
  m_exodusNode3FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Exodus_Node_3");
  m_exodusNode4FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Exodus_Node_4");
  m_exodusNode5FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Exodus_Node_5");
  m_exodusNode6FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Exodus_Node_6");
  m_exodusNode7FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Exodus_Node_7");
  m_exodusNode8FieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Exodus_Node_8");

  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_errorFieldId);
  m_fieldIds.push_back(m_globalErrorFieldId);
  m_fieldIds.push_back(m_exodusNode1FieldId);
  m_fieldIds.push_back(m_exodusNode2FieldId);
  m_fieldIds.push_back(m_exodusNode3FieldId);
  m_fieldIds.push_back(m_exodusNode4FieldId);
  m_fieldIds.push_back(m_exodusNode5FieldId);
  m_fieldIds.push_back(m_exodusNode6FieldId);
  m_fieldIds.push_back(m_exodusNode7FieldId);
  m_fieldIds.push_back(m_exodusNode8FieldId);
}

//! Destructor.
PeridigmNS::Compute_Error::~Compute_Error(){}

//! Fill the energy vectors
int PeridigmNS::Compute_Error::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks  ) const
{
  double globalError = 0.0;

  for(std::vector<Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){

    int numOwnedPoints = blockIt->numPoints();

    double *modelCoord, *coord, *error;
    double *node1, *node2, *node3, *node4, *node5, *node6, *node7, *node8;
    blockIt->getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoord);
    blockIt->getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coord);
    blockIt->getData(m_errorFieldId, PeridigmField::STEP_NONE)->ExtractView(&error);
    blockIt->getData(m_exodusNode1FieldId, PeridigmField::STEP_NONE)->ExtractView(&node1);
    blockIt->getData(m_exodusNode2FieldId, PeridigmField::STEP_NONE)->ExtractView(&node2);
    blockIt->getData(m_exodusNode3FieldId, PeridigmField::STEP_NONE)->ExtractView(&node3);
    blockIt->getData(m_exodusNode4FieldId, PeridigmField::STEP_NONE)->ExtractView(&node4);
    blockIt->getData(m_exodusNode5FieldId, PeridigmField::STEP_NONE)->ExtractView(&node5);
    blockIt->getData(m_exodusNode6FieldId, PeridigmField::STEP_NONE)->ExtractView(&node6);
    blockIt->getData(m_exodusNode7FieldId, PeridigmField::STEP_NONE)->ExtractView(&node7);
    blockIt->getData(m_exodusNode8FieldId, PeridigmField::STEP_NONE)->ExtractView(&node8);

    for(int localId=0 ; localId<numOwnedPoints ; localId++){

      // Nodes in the original Exodus hex determine the limits of integration
      std::vector<double*> n(8);
      n[0] = &node1[3*localId];
      n[1] = &node2[3*localId];
      n[2] = &node3[3*localId];
      n[3] = &node4[3*localId];
      n[4] = &node5[3*localId];
      n[5] = &node6[3*localId];
      n[6] = &node7[3*localId];
      n[7] = &node8[3*localId];

      // Determine the limits of integration
      double x1(DBL_MAX), x2(-DBL_MAX), y1(DBL_MAX), y2(-DBL_MAX), z1(DBL_MAX), z2(-DBL_MAX);
      double x, y, z;
      for(int i=0 ; i<8 ; ++i){
        x = n[i][0];
        y = n[i][1];
        z = n[i][2];
        if(x < x1)
          x1 = x;
        if(x > x2)
          x2 = x;
        if(y < y1)
          y1 = y;
        if(y > y2)
          y2 = y;
        if(z < z1)
          z1 = z;
        if(z > z2)
          z2 = z;
      }

      // Perform a check to make sure the element is rectangular
      int numX1(0), numX2(0), numY1(0), numY2(0), numZ1(0), numZ2(0);
      double tol = 1.0e-10;
      for(int i=0 ; i<8 ; ++i){
        x = n[i][0];
        y = n[i][1];
        z = n[i][2];
        if(std::abs(x - x1) < tol)
          numX1 += 1;
        if(std::abs(x - x2) < tol)
          numX2 += 1;
        if(std::abs(y - y1) < tol)
          numY1 += 1;
        if(std::abs(y - y2) < tol)
          numY2 += 1;
        if(std::abs(z - z1) < tol)
          numZ1 += 1;
        if(std::abs(z - z2) < tol)
          numZ2 += 1;
      }
      bool isRectangular = (numX1 == 4 && numX2 == 4 && numY1 == 4 && numY2 == 4 && numZ1 == 4 && numZ2 == 4);
      if(!isRectangular){
	std::stringstream ss;
	ss << "**** Error:  Error compute class detected non-rectangular element." << std::endl;
	for(int i=0 ; i<8 ; ++i)
	  ss << "  (" << n[i][0] << ", " << n[i][1] << ", " << n[i][2] << ")" << std::endl;
	ss << "  numX1 = " << numX1 << std::endl;
	ss << "  numX2 = " << numX2 << std::endl;
	ss << "  numY1 = " << numY1 << std::endl;
	ss << "  numY2 = " << numY2 << std::endl;
	ss << "  numZ1 = " << numZ1 << std::endl;
	ss << "  numZ2 = " << numZ2 << std::endl;
	TEUCHOS_TEST_FOR_EXCEPT_MSG(!isRectangular, ss.str());
      }

      // Model coordinates for the node at the center of the element
      double modelCoordX = modelCoord[3*localId];
      double modelCoordY = modelCoord[3*localId+1];
      double modelCoordZ = modelCoord[3*localId+2];

      // Current coordinates for the node at the center of the element
      double coordX = coord[3*localId];
      double coordY = coord[3*localId+1];
      double coordZ = coord[3*localId+2];

      // Displacement
      double dispX = coordX - modelCoordX;
      double dispY = coordY - modelCoordY;
      double dispZ = coordZ - modelCoordZ;

      // Compute the error
      double errorVal = computeError_1(dispX, dispY, dispZ, x1, x2, y1, y2, z1, z2);

      error[localId] = errorVal;
      globalError += errorVal;
    }
  }

  // Sum the global error across processors
  double temp = globalError;
  epetraComm->SumAll(&temp, &globalError, 1);

  // Store global value
  Teuchos::RCP<Epetra_Vector> data = blocks->begin()->getData(m_globalErrorFieldId, PeridigmField::STEP_NONE);
  (*data)[0] = globalError;

  return(0);
}

double PeridigmNS::Compute_Error::computeError_1(double peridynamicSolutionX,
                                                 double peridynamicSolutionY,
                                                 double peridynamicSolutionZ,
                                                 double limitOfIntegration_X_lb,
                                                 double limitOfIntegration_X_ub,
                                                 double limitOfIntegration_Y_lb,
                                                 double limitOfIntegration_Y_ub,
                                                 double limitOfIntegration_Z_lb,
                                                 double limitOfIntegration_Z_ub) const
{
  // The displacement function is:
  //   u_x = U11*x*x
  //   u_y = 0
  //   u_z = 0

  double U11 = 1.0;

  // These are the displacement components as computed by the peridynamic simulation
  double uhx = peridynamicSolutionX;
  double uhy = peridynamicSolutionY;
  double uhz = peridynamicSolutionZ;

  // The limits of integration for a rectangular hex element
  double x1 = limitOfIntegration_X_lb;
  double x2 = limitOfIntegration_X_ub;
  double y1 = limitOfIntegration_Y_lb;
  double y2 = limitOfIntegration_Y_ub;
  double z1 = limitOfIntegration_Z_lb;
  double z2 = limitOfIntegration_Z_ub;

  double h = x2 - x1;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(std::abs(y2 - y1 - h) > 1.0e12, "**** Error:  Error compute class detected non-square element.");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(std::abs(z2 - z1 - h) > 1.0e12, "**** Error:  Error compute class detected non-square element.");

  double x = x1 + 0.5*(x2 - x1);

  // The expression for the displacement error integrated over the volume of the element was derived using Mathematica
  //double error = -1.*(std::pow(uhz,2)*x1 + 0.2*std::pow(U11,2)*std::pow(x1,5) + std::pow(uhx,2)*(x1 - 1.*x2) + std::pow(uhy,2)*(x1 - 1.*x2) - 1.*std::pow(uhz,2)*x2 - 0.2*std::pow(U11,2)*std::pow(x2,5) + U11*uhx*(-0.6666666666666666*std::pow(x1,3) + 0.6666666666666666*std::pow(x2,3)))*(y1 - 1.*y2)*(z1 - 1.*z2);

  //double error = (std::pow(uhx,2)*(-x1 + x2) + (2*U11*uhx*(std::pow(x1,3) - std::pow(x2,3)))/3. + (std::pow(U11,2)*(-std::pow(x1,5) + std::pow(x2,5)))/5.)*(y1 - y2)*(z1 - z2);

  //WRONG  double error = (uhx*uhx + uhy*uhy + uhz*uhz)*h*h*h - (2.0/3.0)*uhx*U11*( pow((x + h/2.0), 3) - pow((x - h/2.0), 3) )*h*h + 0.25*U11*U11*( pow((x + h/2.0), 4) - pow((x - h/2.0), 4) )*h*h;

  double error = (uhx*uhx + uhy*uhy + uhz*uhz)*h*h*h - (2.0/3.0)*uhx*U11*( pow((x + h/2.0), 3) - pow((x - h/2.0), 3) )*h*h + 0.2*U11*U11*( pow((x + h/2.0), 5) - pow((x - h/2.0), 5) )*h*h;

  return error;
}
