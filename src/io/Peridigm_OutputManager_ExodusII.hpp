/*! \file Peridigm_OutputManager_ExodusII.hpp */
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
#ifndef PERIDIGM_OUTPUTMANAGER_EXODUSII_HPP
#define PERIDIGM_OUTPUTMANAGER_EXODUSII_HPP

#include <map>

#include <Peridigm_OutputManager.hpp>

#include <Teuchos_ParameterList.hpp>

// Forward declaration
namespace PeridigmNS {
  class Peridigm; 
} 

namespace PeridigmNS {
  
  class OutputManager_ExodusII: public PeridigmNS::OutputManager {
    
  public:
    
    //! Basic constructor.
    OutputManager_ExodusII(const Teuchos::RCP<Teuchos::ParameterList>& params, 
                          PeridigmNS::Peridigm *peridigm_,
                          Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks);
    
    //! Destructor.
    virtual ~OutputManager_ExodusII();

    //! Write data to disk
    virtual void write(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks, double);

  private:
    
    //! Copy constructor.
    OutputManager_ExodusII( const OutputManager& OM );
    
    //! Assignment operator.
    OutputManager_ExodusII& operator=( const OutputManager& OM );

    //! Initialize a new exodus database
    void initializeExodusDatabase(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks);

    //! Initialize a new exodus database that contains only global data
    void initializeExodusDatabaseWithOnlyGlobalData(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks);

    //! Error & Warning reporting tool for calls to ExodusII API
    void reportExodusError(int errorCode, const char *methodName, const char *exodusMethodName);

    //! Valid Teuchos::ParameterList 
    Teuchos::ParameterList getValidParameterList();

    //! Parent pointer
    PeridigmNS::Peridigm *peridigm;

    // Filename of current exodus database
    std::ostringstream filename;

    //! Exodus file handle
    int file_handle;

    //! Index of number of timesteps data actually written to exodus file
    int exodusCount;

    //! Index of first plot dump step to Exodus file
    int firstOutputStep;
    
    //! Index of last plot dump step to Exodus file
    int lastOutputStep;

    //! Flag indicating if this is the first call to initializeExodusDatabase
    bool initializeExodusDatabaseCalled;

    //! Word sizes for IO and CPU
    int CPU_word_size, IO_word_size;

    //! Flag indicating Exodus databases that contain only global data
    bool globalDataOnly;

    //! Map from global output field name to integer. Exodus uses an integer (1..k)  to index the output fields
    std::map <string, int> global_output_field_map;

    //! Map from nodal output field name to integer. Exodus uses an integer (1..k)  to index the output fields
    std::map <string, int> node_output_field_map;

    //! Map from element field name to integer. Exodus uses an integer (1..k)  to index the output fields
    std::map <string, int> element_output_field_map;

    //! Field id for processor id.
    int procNumFieldId;

    //! Field id for element id.
    int elementIdFieldId;
  };
  
}
 
#endif //PERIDIGM_OUTPUTMANAGER_EXODUSII_HPP
