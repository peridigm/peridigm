/*! \file Peridigm_Timer.cpp */

#include "Peridigm_Timer.hpp"
#include <iostream>

using namespace std;

PeridigmNS::Timer& PeridigmNS::Timer::self() {
  static Timer timer;
  return timer;
}

void PeridigmNS::Timer::printTimingData(ostream &out){
  
  //! \todo Parallel communication to allow for printing of min, max, and total values; currently just printing stats from calling processor.

  out << "Timing data (seconds):" << endl;

  for(map<string, TimeKeeper>::const_iterator it=timers.begin() ; it!=timers.end() ; it++){
    const string& name = it->first;
    const TimeKeeper& timeKeeper = it->second;
    out << "  " << name << ":  " << timeKeeper.getElapsedTime() << endl;
  }

}
