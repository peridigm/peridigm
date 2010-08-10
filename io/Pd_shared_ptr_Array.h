/*
 * Pd_shared_ptr.h
 *
 *  Created on: Mar 30, 2010
 *      Author: jamitch
 */

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
