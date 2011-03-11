/*
 * Array.h
 *
 *  Created on: Mar 11, 2011
 *      Author: jamitch
 */

#ifndef ARRAY_H_
#define ARRAY_H_


namespace UTILITIES {

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
	Array(size_t length) : aPtr(new T[length], Deleter()), size(length) {}
	size_t get_size() const { return size; }
	T* get() { return aPtr.get(); }
	const T* get() const { return aPtr.get(); }

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
	size_t size;

};

}


#endif /* ARRAY_H_ */
