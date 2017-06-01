/*! \file Peridigm_UserDefinedRKPMKernel.cpp */

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

#include "Peridigm_UserDefinedRKPMKernel.hpp"
#include "Peridigm_Field.hpp"
#include "rkpmshapefunction.h"

using namespace std;

PeridigmNS::UserDefinedRKPMKernel::UserDefinedRKPMKernel(const Teuchos::ParameterList& params)
  : RKPMKernel(params), 
    m_modelCoordinatesFieldId(-1), 
    m_coordinatesFieldId(-1),  
    m_bondDamageFieldId(-1), 
    m_OwnRKPMShapeFunctionFieldId(-1), 
    m_BondRKPMShapeFunctionFieldId(-1)
{
  functionkernel = params.get<string>("User Defined Kernel");
  m_RKPMSupport = params.get<double>("RKPM Support");
  m_RKPMBasisOrder = params.get<int>("RKPM Basis Order");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Coordinates");
  m_bondDamageFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Bond_Damage");
  m_OwnRKPMShapeFunctionFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Own_RKPM_Shape_Function");
  m_BondRKPMShapeFunctionFieldId   = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Bond_RKPM_Shape_Function");
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_OwnRKPMShapeFunctionFieldId);
  m_fieldIds.push_back(m_BondRKPMShapeFunctionFieldId);

}

PG_RuntimeCompiler::Function PeridigmNS::UserDefinedRKPMKernel::rtcFunction(2, "rtcUserDefinedRKPMKernel");

int PeridigmNS::UserDefinedRKPMKernel::CountKernelDefined = 0;

PeridigmNS::UserDefinedRKPMKernel::~UserDefinedRKPMKernel()
{
}


void PeridigmNS::UserDefinedRKPMKernel::evaluateParserKernel()
{
    if (CountKernelDefined < 1)
    {
    
      rtcFunction.addVar("double", "zeta");
      rtcFunction.addVar("double", "phi");

      string rtcFunctionString = functionkernel;
      if(rtcFunctionString.find("phi") == string::npos)
        rtcFunctionString = "phi = " + rtcFunctionString;

      bool success = rtcFunction.addBody(rtcFunctionString);
      if (!success)
      {
        string msg = "\n**** error:  rtcFunction.addbody(functiondkernel) returned nonzero error code in UserDefinedRKPMKernel::evaluateparserkernel().\n";
        msg += "**** " + rtcFunction.getErrors() + "\n";
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!success, msg);
      }
      success = rtcFunction.varValueFill(1, 0.0);

      CountKernelDefined++;
    }
}

double PeridigmNS::UserDefinedRKPMKernel::userdefinedkernel(const double zeta)
{
    double phi(0.0);
    bool success = rtcFunction.varValueFill(0, zeta);
    if(success)
      success = rtcFunction.execute();
    else
    { 
        string msg = "\n*** Error in UserDefinedRKPMKernel::FillZeta().\n";
        msg += "****" + rtcFunction.getErrors() + "\n";
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!success, msg);
    }
    if(success)
    {
      phi = rtcFunction.getValueOfVar("phi");
    }
    else
    {
        string msg = "\n*** Error in UserDefinedRKPMKernel::CalculatePhi().\n";
        msg += "****" + rtcFunction.getErrors() + "\n";
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!success, msg);
    }
    return phi;
}

void
PeridigmNS::UserDefinedRKPMKernel::initialize(const double dt,
                                              const int numOwnedPoints,
                                              const int* ownedIDs,
                                              const int* neighborhoodList,
                                              PeridigmNS::DataManager& dataManager) const
{
  double *xOverlap, *yOverlap, *bondDamage, *OwnRKPMShapeFunctionValues, *BondRKPMShapeFunctionValues;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
  //dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&yOverlap);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_OwnRKPMShapeFunctionFieldId, PeridigmField::STEP_N)->ExtractView(&OwnRKPMShapeFunctionValues);
  dataManager.getData(m_BondRKPMShapeFunctionFieldId, PeridigmField::STEP_N)->ExtractView(&BondRKPMShapeFunctionValues);
  RKPM_EVALUATION::computeShapeFunction(xOverlap,bondDamage,numOwnedPoints,neighborhoodList,OwnRKPMShapeFunctionValues,BondRKPMShapeFunctionValues,m_RKPMSupport,m_RKPMBasisOrder,userdefinedkernel);

}

void
PeridigmNS::UserDefinedRKPMKernel::computeRKPMShapeFunction(const double dt,
                                                            const int numOwnedPoints,
                                                            const int* ownedIDs,
                                                            const int* neighborhoodList,
                                                            PeridigmNS::DataManager& dataManager) const
{
  double *xOverlap, *yOverlap, *bondDamage, *OwnRKPMShapeFunctionValues, *BondRKPMShapeFunctionValues;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
  //dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&yOverlap);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_OwnRKPMShapeFunctionFieldId, PeridigmField::STEP_NP1)->ExtractView(&OwnRKPMShapeFunctionValues);
  dataManager.getData(m_BondRKPMShapeFunctionFieldId, PeridigmField::STEP_NP1)->ExtractView(&BondRKPMShapeFunctionValues);
  RKPM_EVALUATION::computeShapeFunction(xOverlap,bondDamage,numOwnedPoints,neighborhoodList,OwnRKPMShapeFunctionValues,BondRKPMShapeFunctionValues,m_RKPMSupport,m_RKPMBasisOrder,userdefinedkernel);
}

void 
PeridigmNS::UserDefinedRKPMKernel::applyRKPMShapeFunction(const double dt,
                                                    const int numOwnedPoints,
                                                    const int* ownedIDs,
                                                    const int* neighborhoodList,
                                                    PeridigmNS::DataManager& dataManager) const
{
  double *xOverlap, *yOverlap, *bondDamage, *OwnRKPMShapeFunctionValues, *BondRKPMShapeFunctionValues;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&yOverlap);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_OwnRKPMShapeFunctionFieldId, PeridigmField::STEP_NP1)->ExtractView(&OwnRKPMShapeFunctionValues);
  dataManager.getData(m_BondRKPMShapeFunctionFieldId, PeridigmField::STEP_NP1)->ExtractView(&BondRKPMShapeFunctionValues);
  RKPM_EVALUATION::applyShapeFunction(xOverlap,yOverlap,numOwnedPoints,neighborhoodList,OwnRKPMShapeFunctionValues,BondRKPMShapeFunctionValues);
}

void 
PeridigmNS::UserDefinedRKPMKernel::updateRKPMShapeFunction(const double dt,
                                                     const int numOwnedPoints,
                                                     const int* ownedIDs,
                                                     const int* neighborhoodList,
                                                     PeridigmNS::DataManager& dataManager) const 
{
  *(dataManager.getData(m_OwnRKPMShapeFunctionFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_OwnRKPMShapeFunctionFieldId, PeridigmField::STEP_N));
  *(dataManager.getData(m_BondRKPMShapeFunctionFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_BondRKPMShapeFunctionFieldId, PeridigmField::STEP_N));
}

