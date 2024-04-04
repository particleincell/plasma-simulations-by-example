/*Field is a container for mesh node data division by volume*/
#ifndef _FIELD_H
#define _FIELD_H

#include <ostream>
#include "Vec.h"
#include "VolumeMesh.h"

template <typename T>
class Field_
{
public:
	/*constructor*/
	Field_(VolumeMesh &vm) : ni{(int)vm.nodes.size()}, vm{vm}
	{
		//allocate memory for a 1D array
		data = new T[ni];
		clear();
	}

	//copy constructor
	Field_(const Field_ &other):
	Field_{other.vm} {
		for (int i=0;i<ni;i++)
			data[i] = other(i);
	}
	
	//move constructor
	Field_(Field_ &&other):
		ni{other.ni}, vm{other.vm} {
			data = other.data;	//steal the data
			other.data = nullptr;	//invalidate
	}

	//move assignment operator
	Field_& operator = (Field_ &&f) {return *this;}

	//destructor: release memory
	~Field_() {
		//don't do anything if data is not allocated (or was moved away)
		if (data==nullptr) return;
		delete[] data;
	}

	//overloaded operator [] to allow direct access to data
	T& operator[] (int i) {return data[i];}

	/*returns data[i][j][k] marked as const to signal no data change*/
	T operator() (int i) const {return data[i];}

	/*sets all values to some scalar*/
	void operator =(double s) {
		for (int i=0;i<ni;i++)
			data[i] = s;
	  }

	/*performs element by element division by another field*/
	void operator /= (const Field_ &other) {
		for (int i=0;i<ni;i++)
		{
			if (other.data[i]!=0)
			  data[i] /= other.data[i];
			else
			  data[i] = 0;
		}
	}

	/*increments values by data from another field*/
	Field_& operator += (const Field_ &other) {
		for (int i=0;i<ni;i++)
			data[i]+=other(i);
		return (*this);
	}

	/*performs element by element division by another field*/
	Field_& operator *= (double s) {
		for (int i=0;i<ni;i++)
			data[i]*=s;
		return (*this);
	}

	//multiplication operator, returns f*s
	friend Field_<T> operator*(double s, const Field_<T>&f) {
		Field_<T> r(f);
		return r*=s;
	}

	//scatters value s at logical coordinate lc
	void scatter(const LCord &lc, T s) {
		Tet &tet = vm.tets[lc.cell_id];
		for (int v=0;v<4;v++) data[tet.con[v]] += lc.L[v]*s;
	};

	//gathers field at lc
	T gather(const LCord &lc) {
		T res;
		Tet &tet = vm.tets[lc.cell_id];
		for (int v=0;v<4;v++) res+=lc.L[v]*data[tet.con[v]];
		return res;
	}

	//incorporates new instantaneous values into a running average
	void updateAverage(const Field_ &I) {
		for (int i=0;i<ni;i++)
					data[i] = (I(i)+ave_samples*data[i])/(ave_samples+1);
		++ave_samples;	//increment number of samples
	}


	/*sets all data to zero*/
	void clear() {(*this)=0;}
	
	template<typename S>
	friend std::ostream& operator<<(std::ostream &out, Field_<S> &f);

protected:
	int ni;		//allocated dimensions
	T *data;	//data held by this field
	VolumeMesh &vm; //volume mesh reference
	int ave_samples = 0;	//number of samples used for averaging
};

/*writes out data to a file stream*/
template<typename T>
std::ostream& operator<<(std::ostream &out, Field_<T> &f)
{
	for (int i=0;i<f.ni;i++)
		out<<f.data[i]<<" ";
	return out;
}

//some typedefs
using Field = Field_<double>;
using FieldI = Field_<int>;
using Field3 = Field_<double3>;
using dvector = std::vector<double>;
using dvector3 = std::vector<double3>;
using ivector  = std::vector<int>;

#endif
