/*! \file Peridigm_Timer.hpp */

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

class Timer {

public:

  ~Timer(){}

  //! Singleton.
  static Timer& self();

  void startTimer(std::string name) { timers[name].start(); }

  void stopTimer(std::string name) { timers[name].stop(); }

  void printTimingData(std::ostream &out);

private:

  //! Private constructor
  Timer(){}

  //! @name Private and unimplemented to prevent use
  //@{
  Timer(const Timer&);
  Timer& operator=(const Timer&);
  //@}  

  //! Timer keeper class.
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

    double getElapsedTime() const { return elapsedTime; }

  private:
    Teuchos::RCP<Epetra_Time> epetraTime;
    double elapsedTime;
};

protected:

  std::map<std::string, TimeKeeper> timers;
};

}

#endif // PERIDIGM_TIMER_HPP
