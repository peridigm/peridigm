/*! \file Peridigm_Timer.cpp */

#include "Peridigm_Timer.hpp"
#include <iostream>
#include <vector>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>

using namespace std;

PeridigmNS::Timer& PeridigmNS::Timer::self() {
  static Timer timer;
  return timer;
}

void PeridigmNS::Timer::printTimingData(ostream &out){

  int count = (int)( timers.size() );
  vector<string> names(count);
  vector<double> times(count);
  vector<double> minTimes(count);
  vector<double> maxTimes(count);
  vector<double> totalTimes(count);
  int i = 0;
  for(map<string, TimeKeeper>::reverse_iterator it=timers.rbegin() ; it!=timers.rend() ; it++){
    names[i] = it->first;
    times[i] = it->second.getElapsedTime();
    i++;
  }

  Teuchos::RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  Teuchos::reduceAll<int, double>(*teuchosComm,Teuchos::REDUCE_MIN,count,&times[0], &minTimes[0]);
  Teuchos::reduceAll<int, double>(*teuchosComm,Teuchos::REDUCE_MAX,count,&times[0], &maxTimes[0]);
  Teuchos::reduceAll<int, double>(*teuchosComm,Teuchos::REDUCE_SUM,count,&times[0], &totalTimes[0]);

  unsigned int nameLength = 0;
  for(unsigned int i=0 ; i<names.size() ; ++i)
    if(names[i].size() > nameLength) nameLength = names[i].size();

  int indent = 17;
  int nProc = teuchosComm->getSize();

  if(nProc > 1 && teuchosComm->getRank() == 0){
    out << "Wallclock Time (seconds):" << endl;
    out << "  ";
    out.width(nameLength + 17); out << "Min";
    out.width(indent); out << right << "Max";
    out.width(indent); out << right << "Ave";
    out << endl;
    out.precision(2);
    for(unsigned int i=0 ; i<names.size() ; ++i){
      out.precision(5);
      out << "  ";
      out.width(nameLength + 2); out << left << names[i];
      out.width(indent); out << right << minTimes[i];
      out.width(indent); out << right << maxTimes[i];
      out.width(indent); out << right << totalTimes[i]/nProc;
      out << endl;
    }
    out << endl;
  }
  else if(nProc == 1){
    out << "Wallclock Time (seconds):" << endl;
    out.precision(2);
    for(unsigned int i=0 ; i<names.size() ; ++i){
      out.precision(5);
      out << "  ";
      out.width(nameLength + 2); out << left << names[i];
      out.width(indent); out << right << minTimes[i];
      out << endl;
    }
    out << endl;
  }
}
