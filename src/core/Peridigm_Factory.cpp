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

void updateIntParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn);
void updateDoubleParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn);
void updateBoolParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn);
void updateStringParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn);
void updateParametersFromTextFile(std::string inputFile, Teuchos::Ptr<Teuchos::ParameterList> peridigmParamsPtr);

PeridigmNS::PeridigmFactory::PeridigmFactory(){}

Teuchos::RCP<PeridigmNS::Peridigm> PeridigmNS::PeridigmFactory::create(const std::string inputFile, const MPI_Comm& peridigmComm)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<Epetra_Comm> comm;
  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(peridigmComm));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif

  // Set application parameters to default values
  Teuchos::RCP<Teuchos::ParameterList> peridigmParams = rcp(new Teuchos::ParameterList());
  setPeridigmParamDefaults(peridigmParams.ptr());
  setSolverParamDefaults(peridigmParams.ptr());

  // Determine if string has xml extension
  bool isXML = false;
  std::string ext;
  if (inputFile.size() > 3) {
    ext = inputFile.substr(inputFile.size() - 3);
    if (ext.compare(string("xml")) == 0)
      isXML = true;
  }

  // Update parameters from file
  Teuchos::Ptr<Teuchos::ParameterList> peridigmParamsPtr(peridigmParams.get());
  if (isXML) {
    // Update parameters with data from xml file
    Teuchos::updateParametersFromXmlFile(inputFile, peridigmParamsPtr);
  }
  else {
    updateParametersFromTextFile(inputFile, peridigmParamsPtr);
  }

  // Create new Peridigm object
  return rcp(new PeridigmNS::Peridigm(comm, peridigmParams));

}

void PeridigmNS::PeridigmFactory::setPeridigmParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> peridigmParams_)
{
  peridigmParams_->set("Verbose", false);
}

void PeridigmNS::PeridigmFactory::setSolverParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> peridigmParams_)
{
  Teuchos::ParameterList& solverParams = peridigmParams_->sublist("Solver");

  // general settings
  solverParams.set("Verbose", false);
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
// function dataType to use, if specified data type else use determinable data type
std::vector<std::string> dataType(std::string lineIn)
{
	std::vector<std::string> vec1;
	boost::match_results<std::string::const_iterator> match;
    boost::regex typeMatch("\".*\"");                                             // matches a possible user defined type
	//boost::regex typeMatch("\"[+-]*[a-zA-Z0-9]+.*[a-zA-Z0-9]*[.]*.*[a-zA-Z0-9]*\"");	// matches a possible user defined type
	if(boost::regex_search(lineIn, match, typeMatch))
	{
        //		std::cout << lineIn << std::endl;
		vec1.push_back(lineIn);
		vec1.push_back("string");
	}else{
		vec1.push_back(lineIn);
		vec1.push_back("default");
	}
	return vec1;
}
// function match_ which seperates each line pased in into name, value and if specified data type
std::vector<std::string> match_(std::string line)
{
	std::string str1;
	std::vector<std::string> name_value;
	boost::match_results<std::string::const_iterator> match;
	boost::regex int_re(" ([-+]*[0-9]+[.]*[a-zA-Z0-9]*[a-zA-Z0-9]*[+-]*[0-9]*)");		// matches integers, double or floats
    boost::regex str_re("\".*\"");                                                // matches strings 
	//boost::regex str_re("\"[+-]*[a-zA-Z0-9]+.*_*\\s*[a-zA-Z0-9]*[.]*_*\\s*.*[a-zA-Z0-9]*\"");		// matches strings
	str1 = boost::regex_replace(dataType(line)[0], int_re, "");
	str1 = boost::regex_replace(str1, str_re, "");
    if(boost::regex_search(line, match, str_re) || boost::regex_search(line, match, int_re)){       
	//if(boost::regex_search(line, match, int_re) || boost::regex_search(line, match, str_re)){
		name_value.push_back(cleanUp(str1));
		name_value.push_back(cleanUp(match.str()));
		name_value.push_back(cleanUp(dataType(line)[1]));
        //		std::cout << match.str() << std::endl;
		return name_value;
	}else{
		name_value.push_back(cleanUp(str1));
		name_value.push_back("");
		name_value.push_back(cleanUp(dataType(line)[1]));
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
	bool results = aprepro.parse_stream(infile,inputFile), tab, spc;
    
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
			for(int i = 0; i < line_with_tabs.size()-1; i++){
				std::cout << line_with_tabs[i] << ", ";
			}
			std::cout << line_with_tabs[line_with_tabs.size() - 1] << std::endl;
			throw ;
		}
	}
	return return_data;
}
// function dataParse is intended to loop through data from preparse and determine listing order and style
std::vector<std::pair<std::string, std::string> > dataParse(std::vector<std::pair< std::string, int > > dataIn)
{
	std::string param, line;
	std::pair<std::string, std::string> dataOut;
	std::vector<std::pair<std::string, std::string> > vecOut;
	for(int i = 1; i < dataIn.size(); i++)
	{
		line = dataIn[i-1].first;
		if((dataIn[i-1].second < dataIn[i].second || match_(dataIn[i-1].first)[1].empty()) && match_(dataIn[i-1].first)[1].empty()){
			param = line;
		}
		dataOut = std::make_pair(line, param);
		vecOut.push_back(dataOut);
	}
	dataOut = std::make_pair(dataIn[dataIn.size()-1].first, param);
	vecOut.push_back(dataOut);
	return vecOut;
}
// function dataReturn returns a vector of pairs which contain a vector filled with the line in the first place and the parameter which contains them
std::vector<std::pair<std::vector<std::string>,std::string> > dataReturn(std::vector<std::pair<std::string, std::string> > parsedData, std::vector<std::pair< std::string, int > > dataIn)
{
	std::string param;
	std::vector<std::string> PLVec;
	std::vector<std::pair<std::vector<std::string>,std::string> > returnVec;
	std::pair<std::vector<std::string>,std::string> param_subPLVec;
	for(int i = 0; i < parsedData.size(); i++)
	{
		PLVec.push_back(parsedData[i].first);
		PLVec.push_back(parsedData[i].second);
		if(dataIn[i].second == 0 && match_(dataIn[i].first)[1].empty()){
			param = dataIn[i].first;
		}
		param_subPLVec = std::make_pair(PLVec,param);
		returnVec.push_back(param_subPLVec);
		PLVec.clear();
	}
	return returnVec;
}
// data sorting and distribution sends the value name and type to either one of the
void dataSort(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn, std::string typeIn)
{
	boost::regex reInt("\\s*[-+]*[0-9]+\\s*");
	boost::regex reDouble("\\s*[-+]*[0-9]+.[0-9]+\\s*");
	boost::regex reFloat("\\s*[-+]*[0-9]+[.]*[0-9]*[eE][+-]*[0-9]+\\s*");
	boost::regex reBool("[Ff]alse|[Tt]rue");
    //	std::cout << nameIn << ", " << typeIn << std::endl;
	if(typeIn.compare("default") == 0){
		if(boost::regex_match (valueIn, reInt)){
			updateIntParameter(listPtr, nameIn, valueIn);
		}
		else if(boost::regex_match (valueIn, reDouble) | boost::regex_match (valueIn, reFloat)){
			updateDoubleParameter(listPtr, nameIn, valueIn);
		}
		else if(boost::regex_match (valueIn, reBool)){
			updateBoolParameter(listPtr, nameIn, valueIn);
		}
	}else{
		if(boost::regex_match (valueIn, reBool)){
			updateBoolParameter(listPtr, nameIn, valueIn);
		}else{
			updateStringParameter(listPtr, nameIn, valueIn);
		}
	}
}
void updateIntParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn)
{
	int num;
	listPtr -> set(nameIn, atoi(valueIn.c_str()));
	num = Teuchos::getParameter<int>(*listPtr.get(), nameIn);
}
void updateDoubleParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn)
{
	double dub;
	listPtr-> set(nameIn, std::atof(valueIn.c_str()));
	dub = Teuchos::getParameter<double>(*listPtr.get(), nameIn);
}
void updateBoolParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn)
{
	bool tf;
	listPtr -> set(nameIn, string_to_bool(valueIn));
	tf = Teuchos::getParameter<bool>(*listPtr.get(), nameIn);
}
void updateStringParameter(Teuchos::Ptr<Teuchos::ParameterList> listPtr, std::string nameIn, std::string valueIn)
{
	std::string str;
	listPtr -> set(nameIn,valueIn);
	str = Teuchos::getParameter<std::string>(*listPtr.get(), nameIn);
}
// function updateParametersFromTextFile function that uses each one of the previous functions to store finally parse data to be sent to Peridigm
void updateParametersFromTextFile(std::string inputFile, Teuchos::Ptr<Teuchos::ParameterList> My_List)
{
	std::string previous_list, previous_sublist;
	std::vector<std::pair<std::vector<std::string>, std::string> > dataIn(dataReturn(dataParse(preParse(inputFile)), preParse(inputFile)));
	for(int h = 0; h < dataIn.size(); h++){
		if(dataIn[h].second.compare(dataIn[h].first[0]) == 0){
			Teuchos::Ptr<Teuchos::ParameterList> paramPtr = ptr(new Teuchos::ParameterList());
			previous_list = dataIn[h].second;
			for(int i = h; i<dataIn.size(); i++){
				if(dataIn[i].first[1].compare(dataIn[i].first[0]) == 0 && isspace(dataIn[i].first[0][0]) && previous_list.compare(dataIn[i].second) == 0 ){
					Teuchos::Ptr<Teuchos::ParameterList> listPtr = ptr(new Teuchos::ParameterList());
					previous_sublist = dataIn[i].first[1];
					for(int j = i; j<dataIn.size(); j++){
						if(dataIn[j].first[1].compare(previous_sublist) == 0 && !match_(dataIn[j].first[0])[1].empty()){
							dataSort(listPtr, match_(dataIn[j].first[0])[0], match_(dataIn[j].first[0])[1], match_(dataIn[j].first[0])[2]);
                            //							std::cout << match_(dataIn[j].first[0])[1] << ", " << match_(dataIn[j].first[0])[2] << std::endl;
						}
					}
					paramPtr -> set( cleanUp(previous_sublist), *listPtr.get());
				}
				else if(dataIn[i].first[1].compare(dataIn[i].second)==0 && previous_list.compare(dataIn[i].second)==0 && !match_(dataIn[i].first[0])[1].empty())
				{
					dataSort(paramPtr, match_(dataIn[i].first[0])[0], match_(dataIn[i].first[0])[1], match_(dataIn[i].first[0])[2]);
				}
			}
			My_List -> set(previous_list, *paramPtr.get());
		}
		else if(!isspace(dataIn[h].first[0][0]) && !match_(dataIn[h].first[0])[1].empty())
		{
			dataSort(My_List, match_(dataIn[h].first[0])[0], match_(dataIn[h].first[0])[1], match_(dataIn[h].first[0])[2]);
		}
	}
	if(debug_capture){
		std::string newFileName = inputFile;
		std::cout << newFileName << std::endl;
		while(newFileName[newFileName.size()-1] != '.'){
			newFileName.erase(newFileName.begin() + newFileName.size() - 1);
		}
		newFileName.append("xml");
		const Teuchos::ParameterList paramlist(*My_List.get());
		const std::string xmlFileName(newFileName);
		Teuchos::writeParameterListToXmlFile( paramlist, xmlFileName);
		std::cout << *My_List.get() << std::endl;
	}
}

