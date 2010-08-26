// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHAL_WORKSET_HPP
#define PHAL_WORKSET_HPP

#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <vector>
#include <map>
#include "Peridigm_NeighborhoodData.hpp"
#include "Peridigm_Material.hpp"
#include "Peridigm_ContactModel.hpp"

namespace PHAL {

/*!
 * \brief A workset defines a set of data for processing, for example in a model 
 *        evaluation; for improved performance, the workset may be specified such
 *        that it fits in memory. 
 */
struct Workset {
  
  Workset() {}

  Teuchos::RCP<const Epetra_Vector> xOverlap;
  Teuchos::RCP<const Epetra_Vector> uOverlap;
  Teuchos::RCP<const Epetra_Vector> vOverlap;
  Teuchos::RCP<Epetra_Vector> forceOverlap;
  Teuchos::RCP<Epetra_Vector> contactForceOverlap;
  Teuchos::RCP<const double> timeStep;
  Teuchos::RCP<const Epetra_Vector> cellVolumeOverlap;
  Teuchos::RCP<const PeridigmNS::NeighborhoodData> neighborhoodData;
  Teuchos::RCP<double> bondData;
  Teuchos::RCP<Epetra_MultiVector> scalarConstitutiveDataOverlap;
  Teuchos::RCP<Epetra_MultiVector> vectorConstitutiveDataOverlap;
  Teuchos::RCP<Epetra_MultiVector> bondConstitutiveData;

  Teuchos::RCP<const PeridigmNS::NeighborhoodData> contactNeighborhoodData;

  // The evaluators need access to the material models
  // For now, we're using a vector of materials
  //! \todo Use Teuchos::ArrayRCP to store materials?
  std::vector< Teuchos::RCP<PeridigmNS::Material> > materials;

  std::vector< Teuchos::RCP<PeridigmNS::ContactModel> > contactModels;

  // MPI ID (debugging)
  int myPID;
};

}

#endif
