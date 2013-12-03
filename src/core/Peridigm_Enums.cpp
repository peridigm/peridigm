/*! \file Peridigm_Enums.cpp */

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

#include <Peridigm_Enums.hpp>
#include <map>
#include <boost/algorithm/string.hpp>

static bool string_maps_created = false;

static std::map<Set_Definition,std::string> set_definition_string;
static std::map<std::string,Set_Definition> string_set_definition;
static std::map<Tensor_Order,std::string> tensor_order_string;
static std::map<std::string,Tensor_Order> string_tensor_order;
static std::map<Spatial_Coordinate,std::string> spatial_coordinate_string;
static std::map<std::string,Spatial_Coordinate> string_spatial_coordinate;
static std::map<Boundary_Condition_Type,std::string> boundary_condition_string;
static std::map<std::string,Boundary_Condition_Type> string_boundary_condition;

void create_string_maps();
void create_string_maps()
{
  if (string_maps_created)
  {
    return;
  }
  string_maps_created = true;

  set_definition_string[ALL_SETS]                                                         = "ALL_SETS";
  set_definition_string[FULL_DOMAIN]                                                      = "FULL_DOMAIN";
  set_definition_string[NO_SUCH_SET_DEFINITION]                                           = "NO_SUCH_SET_DEFINITION";

  tensor_order_string[SCALAR]                                                             = "SCALAR";
  tensor_order_string[VECTOR]                                                             = "VECTOR";
  tensor_order_string[TENSOR]                                                             = "TENSOR";
  tensor_order_string[NO_SUCH_TENSOR_ORDER]                                               = "NO_SUCH_TENSOR_ORDER";

  spatial_coordinate_string[X]                                                            = "X";
  spatial_coordinate_string[Y]                                                            = "Y";
  spatial_coordinate_string[Z]                                                            = "Z";
  spatial_coordinate_string[NO_SUCH_SPATIAL_COORDINATE]                                   = "NO_SUCH_SPATIAL_COORDINATE";

  boundary_condition_string[PRESCRIBED_DISPLACEMENT]                                      = "PRESCRIBED_DISPLACEMENT";
  boundary_condition_string[PRESCRIBED_TEMPERATURE]                                       = "PRESCRIBED_TEMPERATURE";
  boundary_condition_string[INITIAL_TEMPERATURE]                                          = "INITIAL_TEMPERATURE";
  boundary_condition_string[INITIAL_DISPLACEMENT]                                         = "INITIAL_DISPLACEMENT";
  boundary_condition_string[INITIAL_VELOCITY]                                             = "INITIAL_VELOCITY";
  boundary_condition_string[BODY_FORCE]                                                   = "BODY_FORCE";
  boundary_condition_string[NO_SUCH_BOUNDARY_CONDITION_TYPE]                              = "NO_SUCH_BOUNDARY_CONDITION_TYPE";


  for (std::map<Set_Definition,std::string>::iterator pos = set_definition_string.begin(); pos != set_definition_string.end(); ++pos)
  {
    string_set_definition[pos->second] = pos->first;
  }
  for (std::map<Tensor_Order,std::string>::iterator pos = tensor_order_string.begin(); pos != tensor_order_string.end(); ++pos)
  {
    string_tensor_order[pos->second] = pos->first;
  }
  for (std::map<Spatial_Coordinate,std::string>::iterator pos = spatial_coordinate_string.begin(); pos != spatial_coordinate_string.end(); ++pos)
  {
    string_spatial_coordinate[pos->second] = pos->first;
  }
  for (std::map<Boundary_Condition_Type,std::string>::iterator pos = boundary_condition_string.begin(); pos != boundary_condition_string.end(); ++pos)
  {
    string_boundary_condition[pos->second] = pos->first;
  }
}

std::string to_string(const Set_Definition & set_definition)
{
  create_string_maps();
  std::map<Set_Definition,std::string>::iterator pos=set_definition_string.find(set_definition);
    if (pos == set_definition_string.end())
    {
      std::stringstream oss;
      oss << "ERROR: Unknown set definition: " << set_definition << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,oss.str());
    }
    return pos->second;
}

Set_Definition to_set_definition(const Teuchos::ParameterList & params)
{
  return to_set_definition(params.get<std::string>("Node Set"));
}

Set_Definition to_set_definition(const std::string & str)
{
  create_string_maps();
  // convert to upper case and trim the input string
  std::string upper_str = str;
  tidy_string(upper_str);

  std::map<std::string,Set_Definition,std::string >::iterator pos=string_set_definition.find(upper_str);
  if (pos == string_set_definition.end())
  {
    return NO_SUCH_SET_DEFINITION;
  }
  return pos->second;
}


std::string to_string(const Tensor_Order & tensor_order)
{
  create_string_maps();
  std::map<Tensor_Order,std::string>::iterator pos=tensor_order_string.find(tensor_order);
    if (pos == tensor_order_string.end())
    {
      std::stringstream oss;
      oss << "ERROR: Unknown tensor order: " << tensor_order << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,oss.str());
    }
    return pos->second;
}

Tensor_Order to_tensor_order(const std::string & str)
{
  create_string_maps();
  // convert to upper case and trim the input string
  std::string upper_str = str;
  tidy_string(upper_str);

  std::map<std::string,Tensor_Order,std::string >::iterator pos=string_tensor_order.find(upper_str);
  if (pos == string_tensor_order.end())
  {
    std::stringstream oss;
    oss << "ERROR: Unknown tensor order: " << upper_str << std::endl;
    oss << "Valid options include: " << std::endl;
    for (std::map<std::string,Tensor_Order,std::string >::iterator it = string_tensor_order.begin(); it != string_tensor_order.end(); ++it)
    {
      oss << "  " << it->first << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,oss.str());
  }
  return pos->second;
}

int to_dimension_size(const Tensor_Order & tensor_order)
{
  create_string_maps();
  int size = -1;
  switch (tensor_order)
  {
  case SCALAR : size = 1;
  break;
  case VECTOR : size = 3; // FIXME assumes 3D
  break;
  case TENSOR : size = 9;  // FIXME assumes 3D
  break;
  case NO_SUCH_TENSOR_ORDER :
  break;
  default :
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"ERROR: unknown tensor order " + to_string(tensor_order));
    break;
  }
  return size;
}

std::string to_string(const Spatial_Coordinate & spatial_coordinate)
{
  create_string_maps();
  std::map<Spatial_Coordinate,std::string>::iterator pos=spatial_coordinate_string.find(spatial_coordinate);
    if (pos == spatial_coordinate_string.end())
    {
      std::stringstream oss;
      oss << "ERROR: Unknown spatial coordinate: " << spatial_coordinate << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,oss.str());
    }
    return pos->second;
}

Spatial_Coordinate to_spatial_coordinate(const Teuchos::ParameterList & params)
{
  if(!params.isParameter("Coordinate")) return X;
  return to_spatial_coordinate(params.get<std::string>("Coordinate"));
}

Spatial_Coordinate to_spatial_coordinate(const std::string & str)
{
  create_string_maps();
  // convert to upper case and trim the input string
  std::string upper_str = str;
  tidy_string(upper_str);

  std::map<std::string,Spatial_Coordinate,std::string >::iterator pos=string_spatial_coordinate.find(upper_str);
  if (pos == string_spatial_coordinate.end())
  {
    std::stringstream oss;
    oss << "ERROR: Unknown spatial coordinate: " << upper_str << std::endl;
    oss << "Valid options include: " << std::endl;
    for (std::map<std::string,Spatial_Coordinate,std::string >::iterator it = string_spatial_coordinate.begin(); it != string_spatial_coordinate.end(); ++it)
    {
      oss << "  " << it->first << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,oss.str());
  }
  return pos->second;
}

int to_index(const Spatial_Coordinate & spatial_coordinate)
{
  create_string_maps();
  int index = -1;
  switch (spatial_coordinate)
  {
  case X : index = 0;
  break;
  case Y : index = 1;
  break;
  case Z : index = 2;
  break;
  case NO_SUCH_SPATIAL_COORDINATE :
  break;
  default :
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"ERROR: unknown spatial coordinate " + to_string(spatial_coordinate));
    break;
  }
  return index;
}

std::string to_string(const Boundary_Condition_Type & bc_type)
{
  create_string_maps();
  std::map<Boundary_Condition_Type,std::string>::iterator pos=boundary_condition_string.find(bc_type);
    if (pos == boundary_condition_string.end())
    {
      std::stringstream oss;
      oss << "ERROR: Unknown boundary condition type: " << bc_type << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,oss.str());
    }
    return pos->second;
}

Boundary_Condition_Type to_boundary_condition_type(const Teuchos::ParameterList & params)
{
  return to_boundary_condition_type(params.get<std::string>("Type"));
}

Boundary_Condition_Type to_boundary_condition_type(const std::string & str)
{
  create_string_maps();
  // convert to upper case and trim the input string
  std::string upper_str = str;
  tidy_string(upper_str);

  std::map<std::string,Boundary_Condition_Type,std::string >::iterator pos=string_boundary_condition.find(upper_str);
  if (pos == string_boundary_condition.end())
  {
    std::stringstream oss;
    oss << "ERROR: Unknown boundary condition type: " << upper_str << std::endl;
    oss << "Valid options include: " << std::endl;
    for (std::map<std::string,Boundary_Condition_Type,std::string >::iterator it = string_boundary_condition.begin(); it != string_boundary_condition.end(); ++it)
    {
      oss << "  " << it->first << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,oss.str());
  }
  return pos->second;
}

void tidy_string(std::string & str){
  boost::to_upper(str);
  boost::algorithm::trim(str);
  while(boost::find_first(str," ")){
    boost::replace_first(str," ","_");
  }
}

bool is_initial(const Boundary_Condition_Type & bc_type){
  const std::string bc_str = to_string(bc_type);
  const std::string test_str = "INITIAL";
  return bc_str.find("INITIAL")!=std::string::npos;
}


