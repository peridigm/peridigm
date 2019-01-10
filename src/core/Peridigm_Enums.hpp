#ifndef PERIDIGM_ENUMS_HPP
#define PERIDIGM_ENUMS_HPP

#include <stdlib.h>
#include <string>
#include <iostream>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

enum Set_Definition
{
  ALL_SETS=0,
  FULL_DOMAIN,
  NO_SUCH_SET_DEFINITION
};

std::string to_string(const Set_Definition & set_definition);
Set_Definition to_set_definition(const std::string & str);
Set_Definition to_set_definition(const Teuchos::ParameterList & params);

enum Spatial_Coordinate
{
  X=0,
  Y,
  Z,
  NO_SUCH_SPATIAL_COORDINATE
};

std::string to_string(const Spatial_Coordinate & spatial_coordinate);
Spatial_Coordinate to_spatial_coordinate(const std::string & str);
Spatial_Coordinate to_spatial_coordinate(const Teuchos::ParameterList & params);
int to_index(const Spatial_Coordinate & spatial_coordinate);

enum Tensor_Order
{
  SCALAR=0,
  VECTOR,
  TENSOR,
  NO_SUCH_TENSOR_ORDER
};

std::string to_string(const Tensor_Order & tensor_order);
Tensor_Order to_tensor_order(const std::string & str);
int to_dimension_size(const Tensor_Order & tensor_order);

enum Boundary_Condition_Type
{
  PRESCRIBED_DISPLACEMENT=0,
  INITIAL_DISPLACEMENT,
  INITIAL_VELOCITY,
  INITIAL_TEMPERATURE,
  INITIAL_FLUID_PRESSURE_U,
  PRESCRIBED_FLUID_PRESSURE_U,
  PRESCRIBED_TEMPERATURE,
  BODY_FORCE,
  NO_SUCH_BOUNDARY_CONDITION_TYPE
};

std::string to_string(const Boundary_Condition_Type & bc_type);
Boundary_Condition_Type to_boundary_condition_type(const std::string & str);
Boundary_Condition_Type to_boundary_condition_type(const Teuchos::ParameterList & params);

std::string trim(std::string const & str);
void tidy_string(std::string & str);
bool is_initial(const Boundary_Condition_Type & bc_type);

#endif /* PERIDIGM_ENUMS_HPP */
