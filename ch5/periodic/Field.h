/*Field is a container for mesh node data division by volume*/
#ifndef _FIELD_H
#define _FIELD_H

#include <ostream>

template <typename T>
struct vec3 {
	vec3 (const T u, const T v, const T w) : d{u,v,w} {}
	vec3 (const T a[3]) : d{a[0],a[1],a[2]} {}
	vec3 (): d{0,0,0} {}
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

template <typename T>
class Field_
{
public:
	
	/*constructor*/
	Field_(int ni, int nj, int nk) :
	ni{ni}, nj{nj}, nk{nk}
	{
		//allocate memory for a 3D array
		data = new T**[ni];
		for (int i=0;i<ni;i++)
		{
			data[i] = new T*[nj];
			for (int j=0;j<nj;j++) data[i][j] = new T[nk];
		}		

		clear();
	}

	//another constructor taking an int3
	Field_(int3 nn) : Field_(nn[0],nn[1],nn[2]) {};

	//copy constructor
	Field_(const Field_ &other):
	Field_{other.ni,other.nj,other.nk} {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k] = other(i,j,k);
	}
	
	//move constructor
	Field_(Field_ &&other):
		ni{other.ni},nj{other.nj},nk{other.nk} {
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

		for (int i=0;i<ni;i++)
		{
			for (int j=0;j<nj;j++)
				delete[] data[i][j];

			delete[] data[i];
		}

		delete[] data;
	}

	//overloaded operator [] to allow direct access to data
	T** operator[] (int i) {return data[i];}

	/*returns data[i][j][k] marked as const to signal no data change*/
	T operator() (int i, int j, int k) const {return data[i][j][k];}

	/*sets all values to some scalar*/
	void operator =(double s) {
		for (int i=0;i<ni;i++)
		  for (int j=0;j<nj;j++)
		   for (int k=0;k<nk;k++)
			data[i][j][k] = s;
	  }

	/*performs element by element division by another field*/
	void operator /= (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++) {
					if (other.data[i][j][k]!=0)
					  data[i][j][k] /= other.data[i][j][k];
				else
					  data[i][j][k] = 0;
			  }
	}

	/*increments values by data from another field*/
	Field_& operator += (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k]+=other(i,j,k);
		return (*this);
	}

	/*performs element by element division by another field*/
	Field_& operator *= (double s) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k]*=s;
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
				for (int k=0;k<f1.nk;k++)
					r[i][j][k] = f1(i,j,k)*f2(i,j,k);
		return r;
	}

	//division of a field by a field of doubles
	friend Field_<T> operator/(const Field_<T>&f, const Field_<double>&d) {
		Field_<T> r(f);
		for (int i=0;i<f.ni;i++)
			for (int j=0;j<f.nj;j++)
				for (int k=0;k<f.nk;k++)
				{
					if (d(i,j,k)!=0)	//check for div by zero
						r[i][j][k] = f(i,j,k)/d(i,j,k);
					else
						r[i][j][k] = 0;
				}
		return r;
	}

	/*returns index for node (i,j,k)*/
	int U(int i, int j, int k) {return k*ni*nj+j*ni+i;}

	/*sets all data to zero*/
	void clear() {(*this)=0;}
	
	/* scatters scalar value onto a field at logical coordinate lc*/
	void scatter(double3 lc, T value)
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		int k = (int)lc[2];
		double dk = lc[2]-k;
		
		data[i][j][k] += (T)value*(1-di)*(1-dj)*(1-dk);
		data[i+1][j][k] += (T)value*(di)*(1-dj)*(1-dk);
		data[i+1][j+1][k] += (T)value*(di)*(dj)*(1-dk);
		data[i][j+1][k] += (T)value*(1-di)*(dj)*(1-dk);
		data[i][j][k+1] += (T)value*(1-di)*(1-dj)*(dk);
		data[i+1][j][k+1] += (T)value*(di)*(1-dj)*(dk);
		data[i+1][j+1][k+1] += (T)value*(di)*(dj)*(dk);
		data[i][j+1][k+1] += (T)value*(1-di)*(dj)*(dk);		
	}

	/* gathers field value at logical coordinate lc*/
	T gather(double3 lc)
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		int k = (int)lc[2];
		double dk = lc[2]-k;
					
		/*gather electric field onto particle position*/
		T val = data[i][j][k]*(1-di)*(1-dj)*(1-dk)+
				data[i+1][j][k]*(di)*(1-dj)*(1-dk)+
				data[i+1][j+1][k]*(di)*(dj)*(1-dk)+
				data[i][j+1][k]*(1-di)*(dj)*(1-dk)+
				data[i][j][k+1]*(1-di)*(1-dj)*(dk)+
				data[i+1][j][k+1]*(di)*(1-dj)*(dk)+
				data[i+1][j+1][k+1]*(di)*(dj)*(dk)+
				data[i][j+1][k+1]*(1-di)*(dj)*(dk);
				
		return val;
	}

	//incorporates new instantaneous values into a running average
	void updateAverage(const Field_ &I) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k] = (I(i,j,k)+ave_samples*data[i][j][k])/(ave_samples+1);
		++ave_samples;	//increment number of samples
	}

	template<typename S>
	friend std::ostream& operator<<(std::ostream &out, Field_<S> &f);
	const int ni,nj,nk;	//allocated dimensions

protected:
	T ***data;	/*data held by this field*/
	int ave_samples = 0;	//number of samples used for averaging
};

/*writes out data to a file stream*/
template<typename T>
std::ostream& operator<<(std::ostream &out, Field_<T> &f)
{
	for (int k=0;k<f.nk;k++,out<<"\n")
		for (int j=0;j<f.nj;j++)
			for (int i=0;i<f.ni;i++) out<<f.data[i][j][k]<<" ";
	return out;
}

//some typedefs
using Field = Field_<double>;
using FieldI = Field_<int>;
using Field3 = Field_<double3>;
using dvector = std::vector<double>;

#endif
