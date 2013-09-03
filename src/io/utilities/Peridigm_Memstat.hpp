/*! \file Peridigm_Memstat.hpp */

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

#ifndef PERIDIGM_MEMSTAT_HPP
#define PERIDIGM_MEMSTAT_HPP

#include <string>
#include <map>
#include <Teuchos_RCP.hpp>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

// This is a very simple class that keeps track of memory use at selected
// points in the code that are usually associated with large allocations
// (i.e. allocating the jacobian, or performing the neighborhood search)
// For more sophisticated profiling, the user should use a tool like valgrind.

namespace PeridigmNS {

//! Simple class for storing and printing memory use stats
class Memstat {

public:
  static Memstat * Instance();

  static Teuchos::RCP<const Epetra_Comm > getComm(){return myComm;}
  static void setComm(Teuchos::RCP<const Epetra_Comm> comm){myComm = comm;}


//  //! Constructor
//  Memstat(){}

//  //! Destructor.
//  ~Memstat(){}

  //! Add a memory stat and catagory to the list
  void addStat(const std::string & description);

  //! Print out the stats
  void printStats();

private:

  //! Constructor
  Memstat(){};

  //! @name Private and unimplemented to prevent use
  //@{
  Memstat(const Memstat&);
  Memstat& operator=(const Memstat&);
  //@}  

  //! Map that associates a description with a stat
  std::map<std::string, unsigned int> stats;

  static Memstat * myMemstatPtr;
  static Teuchos::RCP<const Epetra_Comm> myComm;

};

}

#endif // PERIDIGM_MEMSTAT_HPP
