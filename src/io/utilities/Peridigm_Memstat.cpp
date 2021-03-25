/*! \file Peridigm_Memstat.cpp */

#include "Peridigm_Memstat.hpp"
#if defined(__APPLE__)
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif
#include <iostream>
#include <vector>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>

using namespace std;

PeridigmNS::Memstat * PeridigmNS::Memstat::myMemstatPtr = 0;
Teuchos::RCP<const Epetra_Comm > PeridigmNS::Memstat::myComm;


PeridigmNS::Memstat * PeridigmNS::Memstat::Instance()
{
  if (!myMemstatPtr)   // Only allow one instance of class to be generated.
    myMemstatPtr = new Memstat;
  return myMemstatPtr;
}


void PeridigmNS::Memstat::addStat(const std::string & description){

  size_t heap_size = 0;

#if defined(__GNUC__) && defined(__linux__)
  static struct mallinfo minfo;
  minfo = mallinfo();
  heap_size = minfo.uordblks + minfo.hblkhd;

#elif defined(__APPLE__)
  malloc_statistics_t t = {0,0,0,0};
  malloc_zone_statistics(NULL, &t);
  heap_size = t.size_in_use;

#endif

//  static struct mallinfo minfo;
//  minfo = mallinfo();
//  size_t heap_size = (unsigned int) minfo.uordblks + (unsigned int) minfo.hblkhd;

  // if the descriptor already exists, update it if the memory use is higher
  std::map<std::string, unsigned int>::iterator it = stats.find(description);
  if (it != stats.end()){
    if(heap_size > it->second)
      it->second = heap_size;
  }
  // otherwise create a new entry
  else
    stats.insert(std::pair<std::string,unsigned int>(description,heap_size));
}

void PeridigmNS::Memstat::printStats(){

  if(myComm->NumProc()== 1){
    if(myComm->MyPID() == 0){
      cout << "Memory Usage (Heap Alloc MB):\n";
      //cout << "  " << left << setw(30) << "Checkpoint Description" << right << setw(12) << "Heap Alloc" << "\n";
      for(std::map<std::string,unsigned int>::iterator it=stats.begin();it!=stats.end();++it){
        // if the name is long, trim it:
        std::string desc = it->first;
        if(desc.length() > 25) desc.resize(25);
        cout << "  " << left << setw(30) <<  desc << right << setw(12) << (double)it->second / (1024.0 * 1024.0) << "\n";
      }
      cout << "\n";
    }
  }
  else{
    int count = (int)( stats.size() );
    vector<std::string> names(count);
    vector<double> times(count);
    vector<double> minTimes(count);
    vector<double> maxTimes(count);
    vector<double> totalTimes(count);
    int i = 0;
    for(map<std::string, unsigned int>::reverse_iterator it=stats.rbegin() ; it!=stats.rend() ; it++){
      std::string desc = it->first; // truncate the name if its too long
      if(desc.length() > 25) desc.resize(25);
      names[i] = desc;
      times[i] = it->second;
      i++;
    }
    Teuchos::RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
    Teuchos::reduceAll<int, double>(*teuchosComm,Teuchos::REDUCE_MIN,count,&times[0], &minTimes[0]);
    Teuchos::reduceAll<int, double>(*teuchosComm,Teuchos::REDUCE_MAX,count,&times[0], &maxTimes[0]);
    Teuchos::reduceAll<int, double>(*teuchosComm,Teuchos::REDUCE_SUM,count,&times[0], &totalTimes[0]);
    if(myComm->MyPID() == 0){

      cout << "Memory Usage (Heap Alloc MB):\n";

      cout << "  " << left << setw(30) << " " << right << setw(12) << "Min" << right << setw(15) << "Max" << right << setw(15) << "Ave" << endl;
      for(unsigned int i=0 ; i<names.size() ; ++i){
        cout << "  " << left << setw(30) << names[i] << right << setw(12) << minTimes[i] / (1024.0 * 1024.0)
                                                     << right << setw(15) << maxTimes[i]/ (1024.0 * 1024.0)
                                                     << right << setw(15) << (totalTimes[i]/myComm->NumProc()) / (1024.0 * 1024.0)<< "\n";
      }
      cout << "\n";
    }
  }
}
