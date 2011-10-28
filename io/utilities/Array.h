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

#ifndef ARRAY_H_
#define ARRAY_H_

//#include <tr1/memory>
//#include <memory>
#include "utilities/MemoryInclude.h"
#include <stdexcept>
#include <cmath>
#include <string>

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
	      T* end() { return raw_ptr + size; }

	T operator[](size_t i) const {
		if(size<=i){
			std::string message("ERROR\n\tArray::operator[](size_t i) const \'i\' out of range.");
			throw std::domain_error(message);
		}
		return raw_ptr[i];
	  }

	T & operator[](size_t i) {
		if(size<=i){
			std::string message("ERROR\n\tArray::operator[](size_t i) \'i\' out of range.");
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

	T l2_norm() const {
		const T* s = aPtr.get();
		const T* e = s + size;
		T norm = 0, val;
		for(; s!=e; s++){
			val = *s;
			norm += val * val;
		}
		return norm;
	}


private:
	shared_ptr<T> aPtr;
	T* raw_ptr;
	size_t size;

};

}


#endif /* ARRAY_H_ */
