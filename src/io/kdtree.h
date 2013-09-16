/*
 * kdtree.h
 *
 *  Created on: Nov 9, 2012
 *      Author: jamitch
 */

#ifndef KDTREE_H_
#define KDTREE_H_
#include <cstddef>
#include <algorithm>
#include <limits>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include "utilities/MemoryInclude.h"

namespace femanica {

using std::numeric_limits;
using std::vector;


using std::tr1::shared_ptr;

template<class value_type, class ordinal_type>
class array {

public:



	struct deleter;

	friend struct deleter;

	struct deleter{
		void operator()(value_type* d) {
			delete [] d;
		}
	};

public:
	array() : ptr(), raw_ptr(0), size(0) {}

	array(ordinal_type length) : ptr(new value_type[length], deleter()), raw_ptr(ptr.get()), size(length) {}

	ordinal_type get_size() const { return size; }

	value_type* get() { return ptr.get(); }
	const value_type* get() const { return ptr.get(); }
	const value_type* end() const { return raw_ptr + size; }
	      value_type* end() { return raw_ptr + size; }

	value_type operator[](ordinal_type i) const {
		if(size<=i){
			std::string message("ERROR\n\tarray::operator[](ordinal_type i) const \'i\' out of range.");
			throw std::domain_error(message);
		}
		return raw_ptr[i];
	}

	value_type & operator[](ordinal_type i) {
		if(size<=i){
			std::string message("ERROR\n\tarray::operator[](ordinal_type i) \'i\' out of range.");
			throw std::domain_error(message);
		}
		return raw_ptr[i];
	}

private:
	shared_ptr<value_type> ptr;
	value_type* raw_ptr;
	ordinal_type size;

};





/*
 * I really want my 'enum' to implicitly convert to an 'int'
 */
enum component {NONE=-1,X=0,Y=1,Z=2};

template<class value_type>
struct span_axis {
	component axis;
	value_type cut;
};

template<class value_type>
class point {
public:
      inline point() { p[0] = 0; p[1] = 0; p[2] = 0; }
      inline point(value_type value) { p[0] = value; p[1] = value; p[2] = value; }
      inline point(value_type a, value_type b, value_type c) { p[0] = a ; p[1] = b ; p[2] = c; }
      inline point(const value_type q[3]) { p[0] = q[0] ; p[1] = q[1] ; p[2] = q[2]; }
//	inline point& operator=(const point&) = default;
//	inline point(const point&) = default;

	inline value_type squared() const {
		value_type x=p[X],y=p[Y],z=p[Z];
		return x*x+y*y+z*z;
	}

	inline value_type operator[](int c) const { return p[c%3]; }

	inline bool operator<=(const point& r) const {
		return (
            (p[X]<=r[X]) &&
            (p[Y]<=r[Y]) &&
            (p[Z]<=r[Z])
            );

	}

	inline bool operator<(const point& r) const {
		return (
            (p[X]<r[X]) &&
            (p[Y]<r[Y]) &&
            (p[Z]<r[Z])
            );
	}

  friend std::ostream& operator<<(std::ostream &os, const point<value_type> &p)  {
    os << p[X] << ", " << p[Y] << ", " << p[Z] << std::endl;
    return os;
  }


private:
	value_type p[3];
};


template<class value_type>
point<value_type> intersect(const point<value_type>& p, const span_axis<value_type>& c){
	/*
	 * be careful modifying this function
	 */
	int i=c.axis;
	int j=(c.axis+1)%3;
	int k=(c.axis+2)%3;
	value_type q[3];
	q[i]=c.cut;
	q[j]=p[j];
	q[k]=p[k];
	return point<value_type>(q);
}

template<class value_type>
class rectangular_range {
public:
	/*
	 * default constructor gives largest possible range
	 */
	rectangular_range()
	:low(-numeric_limits<value_type>::max()),high(numeric_limits<value_type>::max()){}

	rectangular_range(const point<value_type>& a, const point<value_type>& b)
	: low(a), high(b) {}

	rectangular_range(const point<value_type>& c, value_type r)
	: low(c[X]-r,c[Y]-r,c[Z]-r), high(c[X]+r,c[Y]+r,c[Z]+r) {}

//	rectangular_range& operator=(const rectangular_range&) = default;
//	rectangular_range(const rectangular_range&) = default;


	rectangular_range left(const span_axis<value_type>& c) const {
		point<value_type> b=femanica::intersect(high,c);
		return rectangular_range(low,b);
	}

	rectangular_range right(const span_axis<value_type>& c) const {
		point<value_type> a=femanica::intersect(low,c);
		return rectangular_range(a,high);
	}

	bool contains(const point<value_type>& p) const {
		return ((low<=p) && (p<=high));
	}

	bool contains(const rectangular_range&r) const {
		return ((low<=r.get_low()) && (r.get_high() <= high));
	}

	bool intersects(const rectangular_range&r) const {
		return (!(high<r.get_low()) && !(r.get_high()<low));
	}

	bool intersects(const rectangular_range&r, component axis) const {
		return (!(high[axis]<r.get_low()[axis]) && !(r.get_high()[axis]<low[axis]));
	}

	point<value_type> get_low() const { return low; }
	point<value_type> get_high() const { return high; }

	friend std::ostream& operator<<(std::ostream &os, const rectangular_range<value_type> &r){
		os << r.get_low() << r.get_high();
		return os;
	}

private:
	const point<value_type> low, high;
};

template<class value_type, class ordinal_type>
class point_array {

public:

	struct point_comparator;
	typedef struct point_comparator {
	public:
		point_comparator(component axis, const value_type *point_data): c(axis), d(point_data) {}
		inline bool operator() (ordinal_type a, ordinal_type b) const { return d[3*a+c] < d[3*b+c]; }
	private:
		const component c;
		const value_type * const d;
	} comparator;

	class iterator;
	friend class iterator;
	class iterator : public std::iterator<std::forward_iterator_tag,value_type,ordinal_type> {
	public:
		iterator() : num_points(0), c(), component_map(),  P(0), p(0), start(0), end(0), data(0) {}
		iterator
		(
				ordinal_type num_points,
				const value_type *point_data,
				array<ordinal_type,ordinal_type> axis_map,
				component axis,
				ordinal_type P=0
		)
		:
			num_points(num_points),
			c(axis),
			component_map(axis_map),
			P(P),
			p(component_map.get()+P),
			start(component_map.get()),
			end(component_map.get()+num_points),
			data(point_data)
		{}
//    iterator& operator=(const iterator&) = default;
//    iterator(const iterator&) = default;

		ordinal_type* map_iterator() const { return p; }
		const ordinal_type* map_end()      const { return end; }
		const ordinal_type* map_start()    const { return start; }
		ordinal_type num_points_from_start() const { return P; }
		ordinal_type num_points_to_end() const { return num_points-P; }

		/*
		 * Need to throw and exception for dereferencing iterator if at "end"
		 */
		value_type operator*()  const { const value_type *d =data+c; return *(d+3*(*p)); }
		value_type operator->() const { const value_type *d = data+c; return *(d+3*(*p)); }
		iterator& operator++()        {if (p!=end) {p++;P++;} return *this;}
		iterator& operator++(int)     { return operator++(); }
		iterator& operator--()        {if (p!=start) {p--;P--;} return *this;}
		iterator& operator--(int)     { return operator--(); }

		iterator operator+(const ordinal_type n) const {
			if (P+n<=num_points){
				return iterator(num_points,data,component_map,c,P+n);
			}
			else {
				/*
				 * Need to throw an exception here; perhaps just return end
				 */
				return *this;
			}
		}

		iterator operator-(const ordinal_type n) const {

			if (P-n>=0) {
				return iterator(num_points ,data,component_map,c,P-n);
			}
			else {
				/*
				 * Need to throw an exception here or perhaps just return end
				 */
				return *this;
			}
		}

		iterator& operator+=(ordinal_type n)    { if (P+n<=num_points) {p+=n;P+=n;} return *this; }
		iterator& operator-=(ordinal_type n)    { if (P-n>=0) {p-=n;P-=n;} return *this; }
		bool operator==( const iterator& r) const  { return p==r.p; }
		bool operator!=( const iterator& r) const  { return p!=r.p; }
		bool operator ()  ( const iterator& left,  const iterator& right) const { return *left < *right; }
		bool operator ()  (value_type left, value_type right) const { return left < right; }

	private:
		ordinal_type num_points;
		component c;
		array<ordinal_type,ordinal_type> component_map;
		ordinal_type P;
		ordinal_type *p;
		const ordinal_type *start, *end;
		const value_type *data;
	};

	iterator begin(component c) const {
		return iterator(num_points,raw_ptr,map[c],c);
	}

	iterator end(component c) const {
		return iterator(num_points,raw_ptr,map[c],c,num_points);
	}

	point_array() : num_points(0), raw_ptr(0), map() {}

	point_array(const value_type* point_data, ordinal_type num_points)
	: num_points(num_points),
	  raw_ptr(point_data),
	  map(3)
	{sort();}



	const array<ordinal_type,ordinal_type>& get_map(component c) const { return map[c]; }
	ordinal_type get_num_points() const { return num_points; }
	const value_type* get() const { return raw_ptr; }

	void sort() {
		map[X]=this->sort(X);
		map[Y]=this->sort(Y);
		map[Z]=this->sort(Z);
	}

	/*
	 * operator access to coordinate data through 'original' ids
	 */
	const value_type* operator[](ordinal_type point_id) const { return (raw_ptr+3*point_id); }

	/*
	 * convenience operator returning a 'point'
	 */
	point<value_type> get_point(ordinal_type point_id) const
	{ return point<value_type>(raw_ptr[3*point_id],raw_ptr[3*point_id+1],raw_ptr[3*point_id+2]); }

private:

	array<ordinal_type,ordinal_type> get_identity_array() const {
		/*
		 * identity used for sorting
		 */

		array<ordinal_type,ordinal_type> ids(num_points);
		ordinal_type *ptr = ids.get();
		ordinal_type *end = ptr + num_points;
		for(ordinal_type j=0;ptr != end; ptr++, j++) {*ptr=j;}
		return ids;

	}

	array<ordinal_type,ordinal_type> sort(component c) const {
		array<ordinal_type,ordinal_type> ids=get_identity_array();
		std::sort(ids.get(),ids.end(),comparator(c,raw_ptr));
		return ids;
	}

	ordinal_type num_points;
	const value_type* raw_ptr;
	array< array<ordinal_type,ordinal_type>, int> map;

};

template<class value_type, class ordinal_type>
struct median {
	typedef typename point_array<value_type, ordinal_type>::iterator iterator;
	ordinal_type id,left_tree_size;
	component axis;
	value_type cut;
	/*
	 * built in copy constructor copies all 3 values of p
	 */
	iterator p[3];
};



template<class value_type, class ordinal_type>
struct span_axis<value_type> get_spanning_axis
(
		const typename point_array<value_type, ordinal_type>::iterator p[3],
		ordinal_type tree_size
){
	/*
	 * intrinsic function 'fabs': ? double versus single
	 */
//	{
//		const typename point_array<value_type, ordinal_type>::iterator e=p[0]+tree_size;
//		typename point_array<value_type, ordinal_type>::iterator i=p[0];
//		typename point_array<value_type, ordinal_type>::iterator j=p[1];
//		typename point_array<value_type, ordinal_type>::iterator k=p[2];
//		for(int I=0;i!=e;I++,i++,j++,k++){
//			std::cout << "span_axis: "<< I << "; " << *i << ", " << *j << ", " << *k << std::endl;
//		}
//		std::cout << "x: start, end: " << *(p[0]) << ", " << *(p[0]+tree_size-1) << std::endl;
//		std::cout << "y: start, end: " << *(p[1]) << ", " << *(p[1]+tree_size-1) << std::endl;
//		std::cout << "z: start, end: " << *(p[2]) << ", " << *(p[2]+tree_size-1) << std::endl;
//	}

	struct span_axis<value_type> a;
	component axis[]={X,Y,Z};
	component c=axis[0];
	typename point_array<value_type, ordinal_type>::iterator i=p[c];
	value_type end=*(i+(tree_size-1));
	value_type start=*i;
	value_type span=fabs(end-start);
//	std::cout << "\t\t" << "span_axis:tree_size = " << tree_size << "\n";
//	std::cout << "\t\t" << c << ":" << "span_axis:start = " << start << "\n";
//	std::cout << "\t\t" << c << ";span_axis:end = " << end << "\n";
//	std::cout << "\t\t" << c << ";span_axis:span = " << span << "\n";
//	std::cout << "\t\t\tchoose span_axis:axis = " << c << "\n";
	for(int n=1;n<3;n++){
		typename point_array<value_type, ordinal_type>::iterator i=p[n];
		end=*(p[n]+(tree_size-1));
		start=*(p[n]);
		value_type m=fabs(end-start);
		if(m>span){
			span=m;
			c=axis[n];
		}
//		std::cout << "\t\t" << n << ":" << "span_axis:start = " << start << "\n";
//		std::cout << "\t\t" << n << ";span_axis:end = " << end << "\n";
//		std::cout << "\t\t" << n << ";span_axis:span = " << m << "\n";
//		std::cout << "\t\t\tchoose span_axis:axis = " << c << "\n";
	}
	a.axis=c;
	a.cut=*(p[c])+span/2.0;
	return a;
}

template<class value_type, class ordinal_type>
const typename point_array<value_type, ordinal_type>::iterator
find_median (const typename point_array<value_type, ordinal_type>::iterator& n, ordinal_type tree_size) {
	value_type min = *n;
	typename point_array<value_type, ordinal_type>::iterator end(n+tree_size-1);
	value_type end_value=*end;
	value_type span=end_value-min;
	value_type value=min+span/2.0;

	typename point_array<value_type, ordinal_type>::iterator median(std::lower_bound(n,end,value));
	typename point_array<value_type, ordinal_type>::iterator adjacent_less(median);
	adjacent_less--;

	if(value>*(adjacent_less) || n==median)
		return median;

	while(value<=*(adjacent_less) && (n!=adjacent_less)){
		adjacent_less--;
	}

	return adjacent_less;
}

template<class value_type, class ordinal_type>
struct median<value_type,ordinal_type>
find_median
(
		const typename point_array<value_type, ordinal_type>::iterator p[3],
		ordinal_type tree_size,
		const point_array<value_type, ordinal_type>& points,
		femanica::array<ordinal_type,ordinal_type>& scratch
){

	struct span_axis<value_type> a = get_spanning_axis<value_type,ordinal_type>(p,tree_size);

	/*
	 * iterator for spanning axis and coordinate data
	 */
	const typename point_array<value_type, ordinal_type>::iterator m=find_median<value_type, ordinal_type>(p[a.axis],tree_size);
	const value_type* m_data=points.get()+a.axis;

	/*
	 * number of points to splitting value
	 */
	ordinal_type left_tree_size=(m.map_iterator()-p[a.axis].map_iterator());

	/*
	 * iterators for other two axes
	 */
	const typename point_array<value_type, ordinal_type>::iterator q=p[(a.axis+1)%3];
	const typename point_array<value_type, ordinal_type>::iterator r=p[(a.axis+2)%3];

	{
		/*
		 * split q and copy into 'scratch' array
		 */

		ordinal_type *left=scratch.get(), *right=scratch.get()+left_tree_size;
		for(ordinal_type *i=q.map_iterator();i!=(q.map_iterator()+tree_size);i++){
			if(m_data[3*(*i)]<a.cut){
				*left=*i;
				left++;
			}
			else{
				*right=*i;
				right++;
			}
		}
	}

	{
		/*
		 * split r and copy into 'q' array
		 */
		ordinal_type *left=q.map_iterator();
		ordinal_type *right=left+left_tree_size;
		for(ordinal_type *i=r.map_iterator();i!=(r.map_iterator()+tree_size);i++){
			if(m_data[3*(*i)]<a.cut){
				*left=*i;
				left++;
			}
			else{
				*right=*i;
				right++;
			}
		}
	}

	/*
	 * Copy 'scratch' values back into q but this requires a little swap.
	 * Values in 'scratch' are 'q values.'
	 * Values in 'q' are 'r values.'
	 */

	for(ordinal_type i=0,*j=q.map_iterator(),*k=r.map_iterator();i<tree_size;i++,j++,k++){
		ordinal_type q_value_tmp=scratch[i];
		ordinal_type r_value_tmp=*j;
		*j=q_value_tmp;
		*k=r_value_tmp;
	}

	struct median<value_type,ordinal_type> _median;
	_median.id=*(m.map_iterator());
	_median.left_tree_size=left_tree_size;
	_median.axis=a.axis;
	_median.cut=a.cut;
	_median.p[(a.axis+0)%3]=p[a.axis];
	_median.p[(a.axis+1)%3]=q;
	_median.p[(a.axis+2)%3]=r;
	return _median;
}


template<class value_type, class ordinal_type>
class kdtree {
public:

	static kdtree get_tree(const value_type* point_data, ordinal_type num_points){
		if(0==num_points || 0==point_data){
			return kdtree();
		}
		return kdtree(point_data,num_points);
	}

	void all_neighbors_within_radius(const value_type *center, value_type R, vector<ordinal_type>& neighbors) const {

		all_neighbors_cube(center,R,neighbors);
		/*
		 * now restrict points to sphere
		 * largest inscribed cube of side 'h'
		 */
		value_type h=value_type(2.0)*R/sqrt(value_type(3.0));

		/*
		 * use for distance calculation
		 */
		value_type R2=R*R;

		/*
		 * cube centered 'center'
		 * NOTE: divide side 'h' by 2.0
		 */
		point<value_type> p(center[0],center[1],center[2]);
		rectangular_range<value_type> H(p,h/2.0);

		ordinal_type *s=scratch.get();
		for(ordinal_type i=0;i<static_cast<ordinal_type>(neighbors.size());i++){
			ordinal_type j=neighbors[i];

			point<value_type> q= points.get_point(j);
			if(H.contains(q)) {
				*s=j;s++;
				continue;
			}
			/*
			 * point is between inside and outside cube
			 */
			point<value_type> d(p[0]-q[0],p[1]-q[1],p[2]-q[2]);
			if(d.squared()<=R2){
				*s=j;s++;
			}
		}
		neighbors.clear();
		neighbors.insert(neighbors.begin(),scratch.get(),s);
	}

	void all_neighbors_cube(const value_type *center, value_type h, vector<ordinal_type>& neighbors) const {

		/*
		 * clear any pre-existing neighbors;
		 * set 'size=0' but leaves the existing 'capacity' intact
		 */
		neighbors.clear();

		/*
		 * create search range
		 */
		point<value_type> p(center[0],center[1],center[2]);
		rectangular_range<value_type> R(p,h);

		/*
		 * special cases
		 */
		if(0==points.get_num_points())
			return;
		else if(1==points.get_num_points()){
			point<value_type> p=points.get_point(0);
			if(R.contains(p)) neighbors.push_back(0);
			return;
		}

		/*
		 * create unbounded range for root
		 */
		rectangular_range<value_type> r;

		/*
		 * start search with 'root' node
		 * depth=0
		 */
		size_t depth=0;
		search(root,r,R,depth,neighbors);
	}

private:
	typedef typename point_array<value_type, ordinal_type>::iterator iterator;

	enum {dimension=3};

	struct node {
		node(): id(-1), axis(), cut(0), left(0), right(0) {}
		ordinal_type id;
		component axis;
		value_type cut;
		node *left;
		node *right;
	};

	/*
	 * coordinates
	 */
	point_array<value_type, ordinal_type> points;

	/*
	 * tree
	 */
	array<node,ordinal_type> tree;

	/*
	 * tree nodes iterator used for tree construction
	 */
	struct node* next_node;

	/*
	 * root node of tree
	 */
	struct node *root;

	/*
	 * scratch array
	 */
	mutable array<ordinal_type,ordinal_type> scratch;



	kdtree(): points(), tree(0), next_node(0), root(0), scratch(0) {}

	kdtree(const value_type* point_data, ordinal_type num_points)
	: points(point_data,num_points),
	  /*
	   * NOTE
	   * Every point exists as a leaf and as a node with the following exceptions:
	   * 1) The 'least leaf' does not exist as a node.
	   * 'least leaf' is found by starting at the tree root and traversing the tree by always
	   * taking the left node until the left node is null;  the node at the end of this traverse
	   * is the 'least leaf.'
	   * 2) If a tree is created with a single point, it is sort of a node and a leaf; of course its
	   * not much of a tree in this case but it will probably occur in practice.
	   * Therefore we need 2*num_points-1 nodes.
	   */
	  tree (2*num_points-1),
	  next_node(tree.get()),
	  root(0),
	  scratch(num_points)
	{
		size_t depth=0;

		const iterator start[]=
		{
				points.begin(X),
				points.begin(Y),
				points.begin(Z)
		};

		if(num_points>1)
			root=make_tree(start,num_points,depth);
		else if(1==num_points){
			root=next_node;
			root->id=0;
			root->axis=X;
			root->cut=point_data[0];
			root->left=0;
			root->right=0;
		}
	}

	void search
	(
			const node* n,
			const rectangular_range<value_type>& r,
			const rectangular_range<value_type>& R,
			size_t depth,
			vector<ordinal_type>& neighbors
	) const
	{

		/*
		 * if n is a leaf
		 * then report as candidate neighbor
		 */
		//		std::cout << "search: depth = " << depth << std::endl;
		if(0==n->left && 0==n->right){
			//			std::cout << "\tFound leaf: " << "n->id = " << n->id << std::endl;
			const value_type *d=points[n->id];
			point<value_type> p(*d,*(d+1),*(d+2));
			if(R.contains(p)){
				//				std::cout << "\t\tFOUND point contained in search range R;" << "n->id = " << n->id << std::endl;
				neighbors.push_back(n->id);
			}
			return;
		}
		/*
		 * Left
		 */

      span_axis<value_type> sa={n->axis,n->cut};
		rectangular_range<value_type> r_left=r.left(sa);
		rectangular_range<value_type> r_right=r.right(sa);
		if (R.contains(r_left)){
			/*
			 * R contains all points below 'n'
			 * report all
			 */
			//			std::cout << "\tSearch range \'contains left\' range; n->id " << n->id << std::endl;
			report_node(n->left,neighbors);

		} else if(R.intersects(r_left,n->axis)){
			/*
			 * R intersects set of points contained below 'n'
			 * search left
			 */
			//			std::cout << "\tSearch range \'intersects left\' range; n->id " << n->id << std::endl;
			//			std::cout << "\tr_left\n\t\t" << r_left.get_low() << "\t\t" << r_left.get_high();
			search(n->left,r_left,R,depth+1,neighbors);
		}

		/*
		 * Right
		 */
		if (R.contains(r_right)){
			/*
			 * R contains all points below 'n'
			 * report all
		     */
			//			std::cout << "\tSearch range \'contains right\' range; n->id " << n->id << std::endl;			return;
			report_node(n->right,neighbors);

		} else if(R.intersects(r_right,n->axis)){
			/*
			 * R intersects set of points contained below 'n'
			 * search right
			 */
			//			std::cout << "\tSearch range \'intersects right\' range; n->id " << n->id << std::endl;
			//			std::cout << "\tr_right\n\t\t" << r_right.get_low() << "\t\t" << r_right.get_high();
			search(n->right,r_right,R,depth+1,neighbors);
		}


	}

	void report_tree(vector<ordinal_type>& neighbors) const {
		/*
		 * clear any pre-existing neighbors;
		 * set 'size=0' but leaves the existing 'capacity' intact
		 */
		neighbors.clear();
		report_node(root,neighbors);
	}

	void report_node(const node* n, vector<ordinal_type>& neighbors) const {
		/*
		 * leaf
		 */
		if(0==n->left && 0==n->right) { neighbors.push_back(n->id); return; }
		/*
		 * not a leaf
		 */
		if(0!=n->left) report_node(n->left,neighbors);
		if(0!=n->right) report_node(n->right,neighbors);
		return;
	}



	struct node* make_tree(const iterator n[3], ordinal_type tree_size, size_t depth){

		if(0==tree_size)
			return 0;

		/*
		 * new node is associated with median
		 */
		struct kdtree<value_type,ordinal_type>::node *root=next_node; next_node++;
		struct median<value_type,ordinal_type> m=find_median<value_type,ordinal_type>(n,tree_size,points,scratch);
		const ordinal_type left_tree_size=m.left_tree_size;
		const ordinal_type right_tree_size=tree_size-left_tree_size;
		const iterator n_right[]={m.p[0]+left_tree_size,m.p[1]+left_tree_size,m.p[2]+left_tree_size};
		//		std::cout << "\tdepth = " << depth << std::endl;
		//		std::cout << "\t\ttree_size: " << tree_size << "\n";
		//		std::cout << "\t\tcut axis: " << m.axis << "\n";
		//		std::cout << "\t\tcut value: " << m.cut << "\n";
		//		std::cout << "\t\tleft_tree_size: " << left_tree_size << "\n";
		//		std::cout << "\t\tright_tree_size: " << right_tree_size << "\n";
		root->id=m.id;
		root->axis=m.axis;
		root->cut=m.cut;
		if(1==left_tree_size){
			struct kdtree<value_type,ordinal_type>::node *left=next_node; next_node++;
			left->id=*(n[m.axis].map_iterator());
			left->axis=NONE;
			left->cut=*(n[m.axis]);
			left->left=0;
			left->right=0;
			root->left=left;
			//			std::cout << "\tnode id: " <<root->id << "\n";
			//			std::cout << "\t\t left id: " << (root->left->id) << "\n";
			//			std::cout << "\t\t tree_size: " << tree_size << "\n";
			//			std::cout << "\t\t left_tree_size: " << left_tree_size << "\n";

		} else root->left=make_tree(n,left_tree_size,depth+1);

		if(1==right_tree_size){
			struct kdtree<value_type,ordinal_type>::node *right=next_node; next_node++;
			right->id=*(n_right[m.axis].map_iterator());
			right->axis=NONE;
			right->cut=*(n_right[m.axis]);
			right->left=0;
			right->right=0;
			root->right=right;
			//			std::cout << "\tnode id: " << root->id << "\n";
			//			std::cout << "\t\tright id: " << (root->right->id) << "\n";
			//			std::cout << "\t\ttree_size: " << tree_size << "\n";
			//			std::cout << "\t\tright_tree_size: " << right_tree_size << "\n";
		} else root->right=make_tree(n_right, right_tree_size, depth+1);
		return root;
	}


};


}




#endif /* KDTREE_H_ */
