/*! \file Peridigm_OutputManager.hpp */
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
#ifndef PERIDIGM_OUTPUTMANAGER_HPP
#define PERIDIGM_OUTPUTMANAGER_HPP

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Peridigm_Block.hpp>
#include <Peridigm_DataManager.hpp>
#include <Peridigm_NeighborhoodData.hpp>

namespace PeridigmNS {

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
    virtual void write(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks, double) = 0;

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
    // ASCII or BINARY?
    std::string outputFormat;
    // Write full neighborlist for each particle?
    bool writeNeighborlist;
    // Filename base
    std::string filenameBase;
    // Parameterlist of user-requested data for output
    Teuchos::RCP<Teuchos::ParameterList> outputVariables;

  private:

    //! Copy constructor.
    OutputManager( const OutputManager& OM );

    //! Assignment operator.
    OutputManager& operator=( const OutputManager& OM );

  };

}

#endif //PERIDIGM_OUTPUTMANAGER_HPP
