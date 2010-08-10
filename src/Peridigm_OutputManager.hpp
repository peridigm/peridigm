/*! \file Peridigm_OutputManager.hpp */
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
#ifndef PERIDIGM_OUTPUTMANAGER_HPP
#define PERIDIGM_OUTPUTMANAGER_HPP

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include "Peridigm_NeighborhoodData.hpp"

namespace Peridigm {
  
  class OutputManager {
    
  public:
    
    //! Basic constructor.
    OutputManager(){};
    
    //! Destructor.
    virtual ~OutputManager() {};

    //! Open file
    virtual void open(){};

    //! Close file
    virtual void close(){};

    //! Write data to disk
    virtual void write(Teuchos::RCP<const Epetra_Vector>, 
					   Teuchos::RCP<const Epetra_MultiVector>,
					   Teuchos::RCP<const NeighborhoodData>,
					   Teuchos::RCP<Teuchos::ParameterList>&) = 0;

  protected:

    //! Number of processors and processor ID
    int numProc;
    int myPID;

    // True if this object writes to disk
    bool iWrite;
    // Number of times write() has been called
    int count;
    // Output frequency
    int frequency;
    // Serial or Parallel?
    bool parallelWrite;
    // ASCII or BINARY?
    string outputFormat;
    // Write full neighborlist for each particle?
    bool writeNeighborlist;
    // Filename base
    string filenameBase;
    // Parameterlist of user-requested data for output
    Teuchos::RCP<Teuchos::ParameterList> materialOutputFields;

  private:
    
    //! Copy constructor.
    OutputManager( const OutputManager& OM );
    
    //! Assignment operator.
    OutputManager& operator=( const OutputManager& OM );

  };
  
}
 
#endif //PERIDIGM_OUTPUTMANAGER_HPP
