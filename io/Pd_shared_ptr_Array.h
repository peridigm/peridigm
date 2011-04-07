/*! \file Pd_shared_ptr_Array.h */

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

#ifndef PD_SHARED_PTR_ARRAY_H_
#define PD_SHARED_PTR_ARRAY_H_
#include <cstddef>
#include <tr1/memory>

template<typename T>
class Pd_shared_ptr_Array {


public:
	struct Deleter;
	friend struct Deleter;
	struct Deleter{
		void operator()(T* d) {
			delete [] d;
		}
	};

	Pd_shared_ptr_Array(std::size_t size=0);
	Pd_shared_ptr_Array(std::size_t size, std::tr1::shared_ptr<T> d);
	Pd_shared_ptr_Array           (const Pd_shared_ptr_Array<T> & copyMe);
	Pd_shared_ptr_Array& operator=(const Pd_shared_ptr_Array<T> & copyMe);
	T* get();
	const T* get() const;
	const T* const end() const;
	std::tr1::shared_ptr<T>& get_shared_ptr();
	std::size_t getSize() const;

private:
	std::tr1::shared_ptr<T> data;
	std::size_t size;

};

template<typename T>
inline Pd_shared_ptr_Array<T>::Pd_shared_ptr_Array(std::size_t size):
data(),
size(size)
{
	data = std::tr1::shared_ptr<T>(new T[size],Deleter());
}

template<typename T>
inline Pd_shared_ptr_Array<T>::Pd_shared_ptr_Array(std::size_t size, std::tr1::shared_ptr<T> d):
data(d),
size(size)
{}

template<typename T>
inline Pd_shared_ptr_Array<T>::Pd_shared_ptr_Array(const Pd_shared_ptr_Array<T> & copyMe)
{
	if(this != &copyMe) {
		data = copyMe.data;
		size = copyMe.size;
	}

}

template<typename T>
inline Pd_shared_ptr_Array<T>& Pd_shared_ptr_Array<T>::operator=(const Pd_shared_ptr_Array<T> & copyMe)
{
	if(this == &copyMe) return *this;
	data = copyMe.data;
	size = copyMe.size;
	return *this;
}


template<typename T>
inline T* Pd_shared_ptr_Array<T>::get() { return data.get(); }

template<typename T>
inline const T* Pd_shared_ptr_Array<T>::get() const { return data.get(); }

template<typename T>
inline const T* const Pd_shared_ptr_Array<T>::end() const { return data.get()+size; }

template<typename T>
inline std::tr1::shared_ptr<T>& Pd_shared_ptr_Array<T>::get_shared_ptr(){ return data; }


template<typename T>
inline std::size_t Pd_shared_ptr_Array<T>::getSize() const  { return size; }


#endif /* PD_SHARED_PTR_ARRAY_H_ */
