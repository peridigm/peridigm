/*! \file Peridigm_OutputManager_VTK_XML.hpp */
// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
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
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************
#ifndef PERIDIGM_OUTPUTMANAGER_VTK_XML_HPP
#define PERIDIGM_OUTPUTMANAGER_VTK_XML_HPP

#include <Peridigm_OutputManager.hpp>
#include <PdVTK.h>

#include <Teuchos_ParameterList.hpp>

namespace PeridigmNS {
  
  class OutputManager_VTK_XML: public PeridigmNS::OutputManager {
    
  public:
    
    //! Basic constructor.
    OutputManager_VTK_XML(const Teuchos::RCP<Teuchos::ParameterList>& params);
    
    //! Destructor.
    virtual ~OutputManager_VTK_XML();

    //! Write data to disk
    virtual void write(Teuchos::RCP<const Epetra_Vector>,
		       Teuchos::RCP<const Epetra_Vector>, 
		       Teuchos::RCP<const Epetra_Vector>, 
		       Teuchos::RCP<const Epetra_Vector>, 
		       Teuchos::RCP<const Epetra_Vector>, 
		       Teuchos::RCP<const Epetra_MultiVector>, 
		       Teuchos::RCP<const NeighborhoodData>, 
		       Teuchos::RCP<Teuchos::ParameterList>&);

  private:
    
    //! Copy constructor.
    OutputManager_VTK_XML( const OutputManager& OM );
    
    //! Assignment operator.
    OutputManager_VTK_XML& operator=( const OutputManager& OM );

    Teuchos::RCP<PdVTK::CollectionWriter> vtkWriter;

    vtkSmartPointer<vtkUnstructuredGrid> grid;

  };
  
}
 
#endif //PERIDIGM_OUTPUTMANAGER_VTK_XML_HPP
