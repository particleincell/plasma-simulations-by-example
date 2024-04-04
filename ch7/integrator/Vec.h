#ifndef _VEC_H
#define _VEC_H

#include <math.h>
#include <fstream>

/*define constants*/
namespace Const
{
	const double EPS_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
	const double QE = 1.602176565e-19;		// C, electron charge
	const double AMU = 1.660538921e-27;		// kg, atomic mass unit
	const double ME = 9.10938215e-31;		// kg, electron mass
	const double K = 1.380648e-23;			// J/K, Boltzmann constant
	const double PI = 3.141592653;			// pi
	const double EvToK = QE/K;				// 1eV in K ~ 11604
}

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
	friend vec3<T> unit(const vec3<T> &v) {if (mag(v)>0) return vec3(v)/mag(v); else return {0,0,0};}

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

#endif
