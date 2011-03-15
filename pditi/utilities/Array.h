/*
 * Array.h
 *
 *  Created on: Mar 11, 2011
 *      Author: jamitch
 */

#ifndef ARRAY_H_
#define ARRAY_H_

#include <tr1/memory>
#include <stdexcept>

namespace UTILITIES {

using std::tr1::shared_ptr;

template<class T>
class Array {

public:
	struct Deleter;

	friend struct Deleter;

	struct Deleter{
		void operator()(T* d) {
			delete [] d;
		}
	};

public:
	Array() : aPtr(), raw_ptr(0), size(0) {}

	Array(size_t length) : aPtr(new T[length], Deleter()), raw_ptr(aPtr.get()), size(length) {}

	Array(size_t length, shared_ptr<T> b) : aPtr(b), raw_ptr(aPtr.get()), size(length) {}

	shared_ptr<T> get_shared_ptr() const { return aPtr; }

	size_t get_size() const { return size; }

	T* get() { return aPtr.get(); }

	const T* get() const { return aPtr.get(); }

	const T* end() const { return raw_ptr + size; }

	T operator[](int i) const {
		if(0>i || size<=i){
			std::string message("ERROR\n\tArray::operator[](int i) const \'i\' out of range.");
			throw std::domain_error(message);
		}
		return raw_ptr[i];
	  }

	T & operator[](int i) {
		if(0>i || size<=i){
			std::string message("ERROR\n\tArray::operator[](int i) \'i\' out of range.");
			throw std::domain_error(message);
		}
		return raw_ptr[i];
	  }

	void deep_copy(const Array<T>& rhs) {

		if(&rhs==this)
			return;

		aPtr=shared_ptr<T>(new T[rhs.get_size()], Deleter());
		size=rhs.get_size();
		T* s = aPtr.get();
		const T* r = rhs.get();
		const T* e = s + size;
		for(; s!=e; s++, r++)
			*s = *r;

	}

	void set(T value) {
		T* s = aPtr.get();
		const T* e = s + size;
		for(; s!=e; s++)
			*s = value;
	}

	void scale(T value) {
		T* s = aPtr.get();
		const T* e = s + size;
		for(; s!=e; s++)
			*s *= value;
	}

private:
	shared_ptr<T> aPtr;
	T* raw_ptr;
	size_t size;

};

}


#endif /* ARRAY_H_ */
