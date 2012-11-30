/*! \file Peridigm_Compute_Nearest_Point_Data.cpp */

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

#include <vector>

#include "Peridigm_Compute_Nearest_Point_Data.hpp"
#include "Peridigm_Field.hpp"

using namespace std;

//! Standard constructor.
PeridigmNS::Compute_Nearest_Point_Data::Compute_Nearest_Point_Data(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                                   Teuchos::RCP<const Epetra_Comm> epetraComm_)
  : Compute(params, epetraComm_), m_verbose(false), m_elementIdFieldId(-1), m_modelCoordinatesFieldId(-1),
    m_variableFieldId(-1), m_outputFieldId(-1)
{
  m_positionX = params->get<double>("X");
  m_positionY = params->get<double>("Y");
  m_positionZ = params->get<double>("Z");
  m_variable = params->get<string>("Variable");
  m_outputLabel = params->get<string>("Output Label");

  if(params->isParameter("Verbose"))
    m_verbose = params->get<bool>("Verbose");

  FieldManager& fieldManager = FieldManager::self();
  m_elementIdFieldId = fieldManager.getFieldId("Element_Id");
  m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  m_variableFieldId = fieldManager.getFieldId(m_variable);
  m_outputFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, m_outputLabel);
  m_fieldIds.push_back(m_elementIdFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_variableFieldId);
  m_fieldIds.push_back(m_outputFieldId);
}

//! Destructor.
PeridigmNS::Compute_Nearest_Point_Data::~Compute_Nearest_Point_Data(){}

void PeridigmNS::Compute_Nearest_Point_Data::initialize( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) {

  // Find the node that is closest to the specified location
  double minDistanceSquared = DBL_MAX;
  double actualX, actualY, actualZ;

  for(std::vector<Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){

    Teuchos::RCP<NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();

    double *x, *id;
    blockIt->getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
    blockIt->getData(m_elementIdFieldId, PeridigmField::STEP_NONE)->ExtractView(&id);

    double distanceSquared;
    for(int iID=0 ; iID<numOwnedPoints ; ++iID){
      distanceSquared = (x[3*iID] - m_positionX)*(x[3*iID] - m_positionX) 
        + (x[3*iID+1] - m_positionY)*(x[3*iID+1] - m_positionY) 
        + (x[3*iID+2] - m_positionZ)*(x[3*iID+2] - m_positionZ);
      if(distanceSquared < minDistanceSquared){
        m_elementId = id[iID];
        actualX = x[3*iID];
        actualY = x[3*iID+1];
        actualZ = x[3*iID+2];
        minDistanceSquared = distanceSquared;
      }
    }
  }

  // \todo Parallelize.
//   // reduction operation for parallel computations

//   vector<int> localIntValues(2), globalIntValues(2);
//   vector<double> localDoubleValues(4), globalDoubleValues(4);

//   localIntValues[0] = countTop;
//   localIntValues[1] = countBottom;
//   localDoubleValues[0] = displacementTop;
//   localDoubleValues[1] = forceTop;
//   localDoubleValues[2] = displacementBottom;
//   localDoubleValues[3] = forceBottom;

//   epetraComm()->SumAll(&localIntValues[0], &globalIntValues[0], 2);
//   epetraComm()->SumAll(&localDoubleValues[0], &globalDoubleValues[0], 4);

//   countTop = globalIntValues[0];
//   countBottom = globalIntValues[0];
//   displacementTop = globalDoubleValues[0];
//   forceTop = globalDoubleValues[1];
//   displacementBottom = globalDoubleValues[2];
//   forceBottom = globalDoubleValues[3];

  if(m_verbose && epetraComm->MyPID() == 0){
    stringstream ss;
    ss << "Nearest Point Data Compute Class:" << endl;
    ss << "  Requested variable: " << m_variable << endl;
    ss << "  Requested location: " << m_positionX << ", " << m_positionY << ", " << m_positionZ << endl;
    ss << "  Closest Element Id: " << m_elementId << endl;
    ss << "  Closest Element Position: " << actualX << ", " << actualY << ", " << actualZ << endl;
    cout << ss.str() << endl;
  }
}

int PeridigmNS::Compute_Nearest_Point_Data::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const {

  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error:  The Nearest_Point_Data compute class is a work in progress.\n");

  return 0;
}
