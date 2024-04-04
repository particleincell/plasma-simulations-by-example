/*Field is a container for mesh node data division by volume*/
#ifndef _FIELD_H
#define _FIELD_H

#include <ostream>

// *** vec2 ****
template <typename T>
struct vec2 {
	vec2 (const T u, const T v) : d{u,v} {}
	vec2 (const T a[2]) : d{a[0],a[1]} {}
	vec2 (): d{0,0} {}
	T& operator[](int i) {return d[i];}
	T operator()(int i) const {return d[i];}
	vec2<T>& operator=(double s) {d[0]=s;d[1]=s;return (*this);}
	vec2<T>& operator+=(vec2<T> o) {d[0]+=o[0];d[1]+=o[1];return(*this);}
	vec2<T>& operator-=(vec2<T> o) {d[0]-=o[0];d[1]-=o[1];return(*this);}
	vec2<T> operator/(double s) {vec2<T>o; o[0]=d[0]/s;o[1]=d[1]/s;return o;}
	vec2<T> operator/=(double s) {d[0]/=s;d[1]/=s;return (*this);}

	//dot product of two vectors
	friend T dot(const vec2<T> &v1, const vec2<T> &v2) {
		T s=0;	for (int i=0;i<2;i++) s+=v1(i)*v2(i);
		return s;	}

	//vector magnitude
	friend T mag(const vec2<T> &v) {return sqrt(dot(v,v));}

	//unit vector
	friend vec2<T> unit(const vec2<T> &v) {return vec2(v)/mag(v);}

protected:
	T d[2];
};

//vec2-vec2 operations
template<typename T>	//addition of two vec3s
vec2<T> operator+(const vec2<T>& a, const vec2<T>& b) {
	return vec2<T> (a(0)+b(0),a(1)+b(1));	}
template<typename T>	//subtraction of two vec2s
vec2<T> operator-(const vec2<T>& a, const vec2<T>& b) {
	return vec2<T> (a(0)-b(0),a(1)-b(1));	}
template<typename T>	//element-wise multiplication of two vec2s
vec2<T> operator*(const vec2<T>& a, const vec2<T>& b) {
	return vec2<T> (a(0)*b(0),a(1)*b(1));	}
template<typename T>	//element wise division of two vec3s
vec2<T> operator/(const vec2<T>& a, const vec2<T>& b) {
	return vec2<T> (a(0)/b(0),a(1)/b(1));	}

//vec2 - scalar operations
template<typename T>		//scalar multiplication
vec2<T> operator*(const vec2<T> &a, T s) {
	return vec2<T>(a(0)*s, a(1)*s);}
template<typename T>		//scalar multiplication 2
vec2<T> operator*(T s,const vec2<T> &a) {
	return vec2<T>(a(0)*s, a(1)*s);}

//output
template<typename T>	//ostream output
std::ostream& operator<<(std::ostream &out, vec2<T>& v) {
	out<<v[0]<<" "<<v[1]<<" 0";	//paraview does not support 2-component arrays
	return out;
}

using double2 = vec2<double>;
using int2 = vec2<int>;

// *** vec3 ***************
template <typename T>
struct vec3 {
	vec3 (const T u, const T v, const T w) : d{u,v,w} {}
	vec3 (const T a[3]) : d{a[0],a[1],a[2]} {}
	vec3 (): d{0,0,0} {}
	vec3 (const vec2<T> &o) : d{o(0),o(1),0} {}  //vec2 initializer
	T& operator[](int i) {return d[i];}
	T operator()(int i) const {return d[i];}
	vec3<T>& operator=(double s) {d[0]=s;d[1]=s;d[2]=s;return (*this);}
	vec3<T>& operator+=(vec3<T> o) {d[0]+=o[0];d[1]+=o[1];d[2]+=o[2];return(*this);}
	vec3<T>& operator-=(vec3<T> o) {d[0]-=o[0];d[1]-=o[1];d[2]-=o[2];return(*this);}
	vec3<T> operator/(double s) {vec3<T>o; o[0]=d[0]/s;o[1]=d[1]/s;o[2]=d[2]/s;return o;}
	vec3<T> operator/=(double s) {d[0]/=s;d[1]/=s;d[2]/=s;return (*this);}

	//dot product of two vectors
	friend T dot(const vec3<T> &v1, const vec3<T> &v2) {
		T s=0;	for (int i=0;i<3;i++) s+=v1(i)*v2(i);
		return s;	}

	//vector magnitude
	friend T mag(const vec3<T> &v) {return sqrt(dot(v,v));}

	//unit vector
	friend vec3<T> unit(const vec3<T> &v) {return vec3(v)/mag(v);}

	//cross product
	friend vec3<T> cross(const vec3<T> &a, const vec3<T> &b) {
		return {a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0)};
	}

protected:
	T d[3];
};

//vec3-vec3 operations
template<typename T>	//addition of two vec3s
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)+b(0),a(1)+b(1),a(2)+b(2));	}
template<typename T>	//subtraction of two vec3s
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)-b(0),a(1)-b(1),a(2)-b(2));	}
template<typename T>	//element-wise multiplication of two vec3s
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)*b(0),a(1)*b(1),a(2)*b(2));	}
template<typename T>	//element wise division of two vec3s
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)/b(0),a(1)/b(1),a(2)/b(2));	}

//vec3 - scalar operations
template<typename T>		//scalar multiplication
vec3<T> operator*(const vec3<T> &a, T s) {
	return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}
template<typename T>		//scalar multiplication 2
vec3<T> operator*(T s,const vec3<T> &a) {
	return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}

//output
template<typename T>	//ostream output
std::ostream& operator<<(std::ostream &out, vec3<T>& v) {
	out<<v[0]<<" "<<v[1]<<" "<<v[2];
	return out;
}

using double3 = vec3<double>;
using int3 = vec3<int>;



// *** Field_ ************
template <typename T>
class Field_
{
public:
	
	/*constructor*/
	Field_(int ni, int nj) :
	ni{ni}, nj{nj}
	{
		//allocate memory for a 2D array
		data = new T*[ni];
		for (int i=0;i<ni;i++)
			data[i] = new T[nj];

		clear();
	}

	//another constructor taking an int2
	Field_(int2 nn) : Field_(nn[0],nn[1]) {};

	//copy constructor
	Field_(const Field_ &other):
	Field_{other.ni,other.nj} {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				data[i][j] = other(i,j);
	}
	
	//move constructor
	Field_(Field_ &&other):	ni{other.ni},nj{other.nj} {
			data = other.data;	//steal the data
			other.data = nullptr;	//invalidate
	}

	//move assignment operator
	Field_& operator = (Field_ &&f) {data=f.data;
				f.data=nullptr; return *this;}

	//destructor: release memory
	~Field_() {
		//don't do anything if data is not allocated (or was moved away)
		if (data==nullptr) return;
		for (int i=0;i<ni;i++)	delete[] data[i];
		delete[] data;
	}

	//overloaded operator [] to allow direct access to data
	T* operator[] (int i) {return data[i];}

	/*returns data[i][j] marked as const to signal no data change*/
	T operator() (int i, int j) const {return data[i][j];}

	/*sets all values to some scalar*/
	void operator =(double s) {
		for (int i=0;i<ni;i++)
		  for (int j=0;j<nj;j++)
		   	data[i][j] = s;
	  }

	/*performs element by element division by another field*/
	void operator /= (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
			{
				if (other.data[i][j]!=0)
				  data[i][j] /= other.data[i][j];
				else
				  data[i][j] = 0;
			}  
	}

	/*increments values by data from another field*/
	Field_& operator += (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				data[i][j]+=other(i,j);
		return (*this);
	}

	/*performs element by element division by another field*/
	Field_& operator *= (double s) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				data[i][j]*=s;
		return (*this);
	}

	//multiplication operator, returns f*s
	friend Field_<T> operator*(double s, const Field_<T>&f) {
		Field_<T> r(f);
		return r*=s;
	}

	//division of a field by a field of doubles
	friend Field_<T> operator*(const Field_<T>&f1, const Field_<T>&f2) {
		Field_<T> r(f1);
		for (int i=0;i<f1.ni;i++)
			for (int j=0;j<f1.nj;j++)
				r[i][j] = f1(i,j)*f2(i,j);
		return r;
	}

	//division of a field by a field of doubles
	friend Field_<T> operator/(const Field_<T>&f, const Field_<double>&d) {
		Field_<T> r(f);
		for (int i=0;i<f.ni;i++)
			for (int j=0;j<f.nj;j++)
			{
				if (d(i,j)!=0)	//check for div by zero
					r[i][j] = f(i,j)/d(i,j);
				else
					r[i][j] = 0;
			}
		return r;
	}

	/*returns index for node (i,j)*/
	int U(int i, int j) {return j*ni+i;}

	/*sets all data to zero*/
	void clear() {(*this)=0;}
	
	/* scatters scalar value onto a field at logical coordinate lc*/
	void scatter(double2 lc, T value)
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		data[i][j] += (T)value*(1-di)*(1-dj);
		data[i+1][j] += (T)value*(di)*(1-dj);
		data[i+1][j+1] += (T)value*(di)*(dj);
		data[i][j+1] += (T)value*(1-di)*(dj);
	}

	/* gathers field value at logical coordinate lc*/
	T gather(double2 lc)
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		/*gather electric field onto particle position*/
		T val = data[i][j]*(1-di)*(1-dj)+
				data[i+1][j]*(di)*(1-dj)+
				data[i+1][j+1]*(di)*(dj)+
				data[i][j+1]*(1-di)*(dj);
				
		return val;
	}

	//incorporates new instantaneous values into a running average
	void updateAverage(const Field_ &I) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
					data[i][j] = (I(i,j)+ave_samples*data[i][j])/(ave_samples+1);
		++ave_samples;	//increment number of samples
	}

	template<typename S>
	friend std::ostream& operator<<(std::ostream &out, Field_<S> &f);
	const int ni,nj;	//allocated dimensions

protected:
	T **data;	/*data held by this field*/
	int ave_samples = 0;	//number of samples used for averaging
};

/*writes out data to a file stream*/
template<typename T>
std::ostream& operator<<(std::ostream &out, Field_<T> &f)
{
	for (int j=0;j<f.nj;j++)
		for (int i=0;i<f.ni;i++) out<<f.data[i][j]<<" ";
	return out;
}

//some typedefs
using Field = Field_<double>;
using FieldI = Field_<int>;
using Field3 = Field_<double3>;
using dvector = std::vector<double>;

#endif
