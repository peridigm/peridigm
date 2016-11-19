/*! \file Peridigm_Factory.cpp */

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

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Peridigm_Factory.hpp"
#include "Peridigm.hpp"

// headers from added functions ********************************************************************
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <boost/regex.hpp>
#include <fstream>
#include "aprepro.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
// *****************************************************************

#ifdef USE_YAML
#include <Teuchos_YamlParser_decl.hpp>
#include <Teuchos_YamlParameterListCoreHelpers.hpp>
#endif

//#ifdef HAVE_MPI
  //#include <Epetra_MpiComm.h>
//#else
  //#include <Epetra_SerialComm.h>
//#endif

using namespace std;

void updateIntParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn);
void updateDoubleParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn);
void updateBoolParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn);
void updateStringParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn);
void updateParametersFromTextFile(std::string inputFile, Teuchos::Ptr<Teuchos::ParameterList> peridigmParamsPtr);

PeridigmNS::PeridigmFactory::PeridigmFactory(){}

Teuchos::RCP<PeridigmNS::Peridigm> PeridigmNS::PeridigmFactory::create(const std::string inputFile,
                                                                       const MPI_Comm& comm,
                                                                       Teuchos::RCP<Discretization> inputPeridigmDiscretization)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

//   Teuchos::RCP<Epetra_Comm> comm;
//   #ifdef HAVE_MPI
//     comm = rcp(new Epetra_MpiComm(peridigmComm));
//   #else
//     comm = rcp(new Epetra_SerialComm);
//   #endif

  // Set application parameters to default values
  Teuchos::RCP<Teuchos::ParameterList> peridigmParams = rcp(new Teuchos::ParameterList());
  setPeridigmParamDefaults(peridigmParams.ptr());

  // Determine if string has xml or yaml extension
  bool isXML = false;
#ifdef USE_YAML
  bool isYAML = false;
#endif
  std::string ext;
  if (inputFile.size() > 3) {
    ext = inputFile.substr(inputFile.size() - 3);
    if (ext.compare(string("xml")) == 0)
      isXML = true;
#ifdef USE_YAML
    else {
    ext = inputFile.substr(inputFile.size() - 4);
    if (ext.compare(string("yaml")) == 0)
      isYAML = true;
    }
#endif
  }

  // Update parameters from file
  Teuchos::Ptr<Teuchos::ParameterList> peridigmParamsPtr(peridigmParams.get());
  if (isXML) {
    // Update parameters with data from xml file
    Teuchos::updateParametersFromXmlFile(inputFile, peridigmParamsPtr);
  }
#ifdef USE_YAML
  else if (isYAML) {

	SEAMS::Aprepro aprepro;
	std::ifstream infile(inputFile.c_str());
    
	// get results from aprepro's parsing
    aprepro.parse_stream(infile,inputFile);     // TODO: Check return value (bool).

    // Update parameters with data from yaml string
    //Teuchos::updateParametersFromYamlFile(inputFile, peridigmParamsPtr);
    Teuchos::updateParametersFromYamlCString((aprepro.parsing_results().str()).c_str(), peridigmParamsPtr);
  }
#endif
  else {
    updateParametersFromTextFile(inputFile, peridigmParamsPtr);

  int rank = 0;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  // Print deprecation warning message
  if(rank == 0){

    std::cout << "WARNING!! Peridigms text file input is deprecated and will be \n removed in a future version.  You may consider installing Trilinos \n with the optional TPL_ENABLE_yaml-cpp and convert the *.peridigm input \n to *.yaml, for a similar markup.  Otherwise use the XML input format.\n\n" << std::endl;
  }

  }

  // Create new Peridigm object
  return rcp(new PeridigmNS::Peridigm(comm, peridigmParams, inputPeridigmDiscretization));

}

Teuchos::RCP<PeridigmNS::Peridigm> PeridigmNS::PeridigmFactory::create(const std::string inputFile,
                                                                       const MPI_Comm& comm)
{
  Teuchos::RCP<Discretization> nullDiscretization;
  return create(inputFile, comm, nullDiscretization);
}

void PeridigmNS::PeridigmFactory::setPeridigmParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> peridigmParams_)
{
  peridigmParams_->set("Verbose", false);
}

// function to update parameters from text file*********************************************************************************************
// string to bool used in convert a string "true" or "false" to a boolean "0" or "1"
bool debug_capture = false;
// string to bool used in convert a string "true" or "false" to a boolean "0" or "1"
bool string_to_bool(std::string value)
{
	bool sts;
	std::transform(value.begin(), value.end(), value.begin(), ::tolower);
	std::istringstream in(value);
	in >> std::boolalpha >> sts;
	return sts;
}
// function cleanUp to clear all spaces before and after tag as well as quotations
std::string cleanUp(std::string dataIn)
{
    //	std::cout << dataIn << std::endl;
	boost::regex spcMatch("^(\\s*)|(\\s*)$");		// matches spaces on either side of the value
	dataIn.erase(std::remove(dataIn.begin(), dataIn.end(), '\"'), dataIn.end());
	dataIn = boost::regex_replace(dataIn, spcMatch, "");
	return dataIn;
}
// function match_ which seperates each line pased in into name, value and if specified data type
std::vector<std::string> match_(std::string line)
{
	std::string str1;
	std::vector<std::string> name_value;
	boost::match_results<std::string::const_iterator> match;
	boost::regex int_re(" ([-+]*[0-9]+[.]*[a-zA-Z0-9]*[a-zA-Z0-9]*[+-]*[0-9]*)");		// matches integers, double or float
	boost::regex ver_re("[Vv][Ee][Rr][Bb][Oo][Ss][Ee]");//---------edit 8/29/13---------//
	boost::regex str_re("\".*\"");                                                // matches strings 
	//a extra comment . boost::regex str_re("\"[+-]*[a-zA-Z0-9]+.*_*\\s*[a-zA-Z0-9]*[.]*_*\\s*.*[a-zA-Z0-9]*\"");		// matches strings
	str1 = boost::regex_replace(line, int_re, "");
	str1 = boost::regex_replace(str1, str_re, "");
	if(boost::regex_search(line, match, str_re) || boost::regex_search(line, match, int_re)){       
	//if(boost::regex_search(line, match, int_re) || boost::regex_search(line, match, str_re)){
		name_value.push_back(cleanUp(str1));
		name_value.push_back(match.str());
        //		std::cout << match.str() << std::endl;
		return name_value;
	}else{
		if( boost::regex_search(line, match, ver_re)) {//------------edit 8/29/13----------//|-->>
			name_value.push_back(match.str());
			name_value.push_back("true");
		}else {						//--------edit e/29/13-----------//<<--|
			name_value.push_back(cleanUp(str1));
			name_value.push_back("");
		}
		return name_value;
	}
}
void debug_scan(std::string line)
{
	boost::regex bug(".*[Dd][Ee][Bb][Uu][Gg].*");
	if(boost::regex_search(line, bug)){
		debug_capture = true;
	}
}
// function preParse to store data with flags that indicate the parameterList storage
std::vector<std::pair< std::string, int > > preParse(std::string inputFile)
{
	std::string line;
	SEAMS::Aprepro aprepro;
	std::ifstream infile(inputFile.c_str());
	std::pair<std::string, int> line_data;
    
	// get results from aprepro's parsing
	bool tab, spc;
    aprepro.parse_stream(infile,inputFile);     // \todo Check return value (bool).

	// string stream to store parsed data
	std::istringstream parsed_stream(aprepro.parsing_results().str());
	
	// counting integers
	int space_count, line_count = 0, debug;
	std::vector<int> line_with_tabs;
	// vector to store parsing results
	std::vector<std::pair<std::string, int > > return_data;
	
	while(std::getline(parsed_stream, line))
	{
		tab = false; spc = false;
		if(!line.empty())
		{
			line_count++;
			space_count = 0; debug = 0;
			while(isspace(line[space_count])){
				space_count++;
			}
			while(line[debug] == ' ' || line[debug] == '\t'){
				debug++;
				if(line[debug] == '\t'){
					tab = true;
				}else if(line[debug] == ' '){
					spc = true;
				}
			}
			if(isalnum(line[space_count])){
				line_data = std::make_pair(line, space_count);
				return_data.push_back(line_data);
			}else{
				debug_scan(line);
			}
			if(tab && spc){
				line_with_tabs.push_back(line_count);
			}
		}
	}
	if(!line_with_tabs.empty()){
		try{
			throw 1;
		}catch(int e){
			std::cout << "Warning: Error Number: " << e << std::endl;
			std::cout << "The Combination of Spaces and Tabs Can Affect Parser Accuracy Check Line(s): ";
			for(unsigned int i = 0; i < line_with_tabs.size()-1; i++){
				std::cout << line_with_tabs[i] << ", ";
			}
			std::cout << line_with_tabs[line_with_tabs.size() - 1] << std::endl;
			throw ;
		}
	}
	return return_data;
}
// this function will create a string, int pair where the int will tell the parameterlist level
std::vector<std::pair<std::string, int > > lev_by_line(std::vector<std::pair< std::string, int > > dataIn)
{
	std::vector<std::pair<std::string, int > > line_lev;
	
	// a map to store levels rather than spacing
	std::map < int, int> spc_to_lev;

	// determine the highest level
	int level_count = 0; int space_count = 0;
	spc_to_lev[space_count] = level_count;
	std::vector<std::pair< std::string, int > >::iterator it = dataIn.begin();
	while(it!=dataIn.end()) {
//		std::cout << it->second << std::endl;
		if(it->second > space_count) {
			space_count = it->second;
			level_count++;
			spc_to_lev[space_count] = level_count;
		}
		it++;
	}

	std::advance(it,-std::distance(dataIn.begin(), dataIn.end()));
	// print out number of levels
//	std::cout << level_count << std::endl;
	// iterate through and find leveling scheme
	while(it!=dataIn.end()) {
		line_lev.push_back(std::make_pair(it->first,spc_to_lev[it->second]) );
//		std::cout << it->first << ", " << spc_to_lev[it->second] << std::endl;
		it++;
	}
	return line_lev;
}

// class that contains various functions hardcoded to produce upto 6 level deep parameterlist
class TextParser {
	public:
		static void ParseMyTextFile (std::vector<std::pair<std::string, int> > list_in, Teuchos::Ptr<Teuchos::ParameterList> My_List) {
			// determine parameterlist level
			Plist4deep(list_in, My_List);
		}
	private:
		static void dataSort(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn)
		{
			boost::regex reInt("\\s*[-+]*[0-9]+\\s*");
			boost::regex reDouble("\\s*[-+]*[0-9]+.[0-9]+\\s*");
			boost::regex reFloat("\\s*[-+]*[0-9]+[.]*[0-9]*[eE][+-]*[0-9]+\\s*");
			boost::regex reBool("\\s*\"*[Ff]alse\"*\\s*|\\s*\"*[Tt]rue\"*\\s*");
			if(boost::regex_match (valueIn, reInt)){
				updateIntParameter(listPtr, nameIn, cleanUp(valueIn));
				//std::cout << nameIn << ", " << valueIn << " is an int" << std::endl;
			}
			else if(boost::regex_match (valueIn, reDouble) | boost::regex_match (valueIn, reFloat)){
				updateDoubleParameter(listPtr, nameIn, cleanUp(valueIn));
				//std::cout << nameIn << ", " << valueIn << " is a double" << std::endl;
			}
			else if(boost::regex_match (valueIn, reBool)){
				updateBoolParameter(listPtr, nameIn, cleanUp(valueIn));
				//std::cout << nameIn << ", " << valueIn << " is a bool" << std::endl;
			}
			else {
				updateStringParameter(listPtr, nameIn, cleanUp(valueIn));
				//std::cout << nameIn << ", " << valueIn << " is a string" << std::endl;
			}
		}
		static void updateIntParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn)
		{
			listPtr -> set(nameIn, atoi(valueIn.c_str()));
			Teuchos::getParameter<int>(*listPtr.get(), nameIn);
		}
		static void updateDoubleParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn)
		{
			listPtr-> set(nameIn, std::atof(valueIn.c_str()));
			Teuchos::getParameter<double>(*listPtr.get(), nameIn);
		}
		static void updateBoolParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn)
		{
			listPtr -> set(nameIn, string_to_bool(valueIn));
			Teuchos::getParameter<bool>(*listPtr.get(), nameIn);
		}
		static void updateStringParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn)
		{
			listPtr -> set(nameIn,valueIn);
			Teuchos::getParameter<std::string>(*listPtr.get(), nameIn);
		}
		// 4 levels
		static void Plist4deep (std::vector<std::pair<std::string, int> > list_in, Teuchos::Ptr<Teuchos::ParameterList> My_List) {
			// 3 strings to store 3 different level parameterlist names
			std::pair<std::string,int> previous_0list_data, 
				previous_1list_data, 
				previous_2list_data, 
				previous_3list_data, 
				previous_4list_data, 
				previous_5list_data;
			// create an iterator of the vector coming in
			std::vector<std::pair<std::string,int> >::iterator it = list_in.begin();

			// initiate first 0 level parameterlist name
			previous_0list_data = *it;

			// iterate through data vector coming in and build parameter list from input deck
			while(it!=list_in.end() and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
				// test whether or not while loop needs to begin and prevent unnecessary creations of parameterlist ptr
				if(it->second > 0) {
					// parameter list to store first level parameter name/value and list
					Teuchos::Ptr<Teuchos::ParameterList> lev1list = ptr(new Teuchos::ParameterList());
					// loop as long as list level is greater than 0 level
					while(it->second > 0 and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
						// determine 2nd level needs to be started
						if(it->second > 1) {
							// parameter list to strore 2nd level parameter name/value/list
							Teuchos::Ptr<Teuchos::ParameterList> lev2list = ptr(new Teuchos::ParameterList());
							// loop as long as entry's level is greater than 1 level
							while(it->second > 1 and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
								// start 3rd level?
								if(it->second > 2) {
									// parameter list for 3level
									Teuchos::Ptr<Teuchos::ParameterList> lev3list = ptr(new Teuchos::ParameterList());
									// loop while level is greater than 2
									while(it->second > 2 and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
										// start 4th level?
										if(it->second > 3) {
											// parameter list for 4 level
											Teuchos::Ptr<Teuchos::ParameterList> lev4list = ptr(new Teuchos::ParameterList());
											// loop while level is greater than 3
											while(it->second > 3 and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
												// start 5th level
												if(it->second > 4) {
													// parameterlist for 5 level
													Teuchos::Ptr<Teuchos::ParameterList> lev5list = ptr(
															new Teuchos::ParameterList()
															);
													while(it->second > 4 and 
															std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
														// start 6th level
														if(it->second > 5) {
															// parameterlist for 6 level
															Teuchos::Ptr<Teuchos::ParameterList>
																lev6list = ptr(new
																		Teuchos::ParameterList()
																	      );
															while(it->second > 5 and 
																std::distance(list_in.begin(),it) 
																< static_cast<int>(list_in.size())) {
																// set parameters
																dataSort(lev6list,
																	match_(it->first)[0],
																	match_(it->first)[1]
																	);
																it++;
															}
															lev5list->set(match_(previous_5list_data.first)[0],
																	*lev6list.get()
																	);
														}
														if(it->second == 5 and std::distance(
																	list_in.begin(),it) < static_cast<int>(list_in.size())) {
															if(!match_(it->first)[1].empty()){
																// set parameters
																dataSort(lev5list,
																	match_(it->first)[0],
																	match_(it->first)[1]
																	);
															}
															previous_5list_data = *it;
															it++;
														}
													} // while > 4
													lev4list->set(match_(previous_4list_data.first)[0],
															*lev5list.get()
														     );
												}
												if(it->second == 4 and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
													if(!match_(it->first)[1].empty()) {
														// set parameters 
														dataSort(lev4list,
															match_(it->first)[0],
															match_(it->first)[1]
															);
													}
													previous_4list_data = *it;
													it++;
												}
											} // while > 3
											lev3list->set(match_(previous_3list_data.first)[0],
													*lev4list.get()
												     );
										}
										if(it->second == 3 and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
											if(!match_(it->first)[1].empty()) {
												// set and verify parameters
												dataSort(lev3list,
													match_(it->first)[0],
													match_(it->first)[1]
													);
											}
											previous_3list_data = *it;
											it++;
										}
									} // while > 2
									lev2list->set(match_(previous_2list_data.first)[0],
											*lev3list.get()
										     );
								}
								if(it->second == 2 and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
									// if entry does have a name/value combo
									if(!match_(it->first)[1].empty()) {
										// set parameter
										dataSort(lev2list,
											match_(it->first)[0],
											match_(it->first)[1]
											);
									}
									previous_2list_data = *it;
//									std::cout << previous_2list_data.first << std::endl;
									it++;
								}
							} // while > 1
							lev1list->set(match_(previous_1list_data.first)[0],
									*lev2list.get()
								     );
						}
						if(it->second == 1 and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
							// if entry does has both name and value
							if(!match_(it->first)[1].empty()) {
								// set and verify parameters
								dataSort(lev1list,
									match_(it->first)[0],
									match_(it->first)[1]
									);
							}
							previous_1list_data = *it;
							it++;
						}
					} // while > 0
					// set parameter list of the master list
					My_List->set(match_(previous_0list_data.first)[0],
							*lev1list.get()
						     );
				}
				if(it->second == 0 and std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {
					// if entry has both name and value
					if(!match_(it->first)[1].empty()) {
						// set and verify parameters
						dataSort(My_List,
							match_(it->first)[0],
							match_(it->first)[1]
							);
					}
					previous_0list_data = *it;
				}
				if(std::distance(list_in.begin(),it) < static_cast<int>(list_in.size())) {it++;}
			} // loop for top level
		}
};
// function updateParametersFromTextFile function that uses each one of the previous functions to store finally parse data to be sent to Peridigm
void updateParametersFromTextFile(std::string inputFile, Teuchos::Ptr<Teuchos::ParameterList> My_List)
{
	TextParser alpha;
	alpha.ParseMyTextFile(lev_by_line(preParse(inputFile)), My_List);
	if(debug_capture){
		std::string newFileName = inputFile;
		while(newFileName[newFileName.size()-1] != '.'){
			newFileName.erase(newFileName.begin() + newFileName.size() - 1);
		}
		newFileName.append("xml");
		const Teuchos::ParameterList paramlist(*My_List.get());
		const std::string xmlFileName(newFileName);
		Teuchos::writeParameterListToXmlFile( paramlist, xmlFileName);
	}
	//std::cout << *My_List.get() << std::endl;
}

