/*! \file Peridigm_Timer.hpp */

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

#ifndef PERIDIGM_TIMER_HPP
#define PERIDIGM_TIMER_HPP

#include <Teuchos_RCP.hpp>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 1
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Time.h>
#include <ostream>

namespace PeridigmNS {

//! Singleton class for performance monitoring; manages a set of TimeKeeper objects.
class Timer {

public:

  //! Destructor.
  ~Timer(){}

  //! Singleton.
  static Timer& self();

  //! Starts specified timer, creates the timer if it does not exist.
  void startTimer(std::string name) { timers[name].start(); }

  //! Stops specified timer.
  void stopTimer(std::string name) { timers[name].stop(); }

  //! Prints out a table of timing data.
  void printTimingData(std::ostream &out);

private:

  //! Private constructor
  Timer(){}

  //! @name Private and unimplemented to prevent use
  //@{
  Timer(const Timer&);
  Timer& operator=(const Timer&);
  //@}  

  //! Time keeper class, operates like a stopwatch.
  class TimeKeeper {

  public:

    TimeKeeper() : elapsedTime(0.0) {
#ifdef HAVE_MPI
      epetraTime = Teuchos::rcp(new Epetra_Time(Epetra_MpiComm(MPI_COMM_WORLD)));
#else
      epetraTime = Teuchos::rcp(new Epetra_Time(Epetra_SerialComm()));
#endif
    }

    void start() {
      epetraTime->ResetStartTime();
    }

    void stop() {
      elapsedTime += epetraTime->ElapsedTime();
    }

    //! Returns the cummulative elapsed time.
    double getElapsedTime() const { return elapsedTime; }

  private:
    Teuchos::RCP<Epetra_Time> epetraTime;
    double elapsedTime;
  };

protected:

  //! Map that associates a name with a TimeKeeper.
  std::map<std::string, TimeKeeper> timers;
};

}

#endif // PERIDIGM_TIMER_HPP
