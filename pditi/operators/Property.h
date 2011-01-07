/*
 * Property.h
 *
 *  Created on: Mar 27, 2010
 *      Author: awesome
 */

#ifndef PROPERTY_H_
#define PROPERTY_H_
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

namespace PdImp {

template<typename T>
class Property {
private:
	T value;
	std::string name;
protected:
	explicit Property(const T v, const std::string& param_name) : value(v), name(param_name) {}
public:
	const T operator()() { return value; }
	const T getValue () const { return value; }
	const std::string& getName() const { return name; }
	friend std::ostream& operator << (std::ostream& out, const Property<T>& v) { out << v.getName() << "value  = " << v.getValue() << std::endl; return out;}
	static std::string toString(T const& value) {
	    std::stringstream sstr;
	    sstr << value;
	    return sstr.str();
	}

};


template<typename T>
class Boundary : public Property<T> {
public:
	virtual ~Boundary(){}
	explicit Boundary(T v, std::string name) : Property<T>(v,name){}
	virtual bool operator<(T right) const = 0;
	virtual bool operator>(T right) const = 0;
};

template<typename T>
class Inclusive : public Boundary<T>{
public:
	explicit Inclusive(T v) : Boundary<T>(v, "= ") {}
	bool operator<(T right) const { return this->getValue()<=right; }
	bool operator>(T right) const { return this->getValue()>=right; }
};

template<typename T>
class Exclusive : public Boundary<T>{
public:
	explicit Exclusive(T v) : Boundary<T>(v, " ") {}
	bool operator<(T right) const { return this->getValue()<right; }
	bool operator>(T right) const { return this->getValue()>right; }
};

template<typename T>
class PositiveNumber : public Property<T> {
public:
	PositiveNumber(T v, const std::string& param_name) : Property<T>(v, param_name) {
		if(v <= 0) {
			std::string message(param_name);
			message += " requires number in the range < Input value = ";
			message += Property<T>::toString(v);
			throw std::domain_error(message);
		}
	}
};

template<typename T>
class DomainNumber : public Property<T> {
public:
	explicit DomainNumber(const Boundary<T>& min, const Boundary<T>& max, T value, const std::string& param_name): Property<T>(value, param_name)
	{
		if(!isValid(min, max, value)){
			std::string message(param_name);
			message += " \"x\" is restricted to the following range: x >";
			message += min.getName();
			message += Property<T>::toString(min.getValue());
			message += " AND x <";
			message += max.getName();
			message += Property<T>::toString(max.getValue());
			message += "\n\tInput value = " + Property<T>::toString(value);
			throw std::domain_error(message);
		}
	}
	static bool isValid(const Boundary<T>& min, const Boundary<T>& max, T value) { return min < value && max > value; }
};


}
#endif /* PROPERTY_H_ */
