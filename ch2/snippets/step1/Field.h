#ifndef _FIELD_H
#define _FIELD_H

//Field.h
class Field {
public:

	//constructor
	Field (int ni, int nj, int nk) :
		ni{ni}, nj{nj}, nk{nk} {
			
			data = new double**[ni];	//allocate ni pointers-to-pointers
			for (int i=0;i<ni;i++)
			{
				data[i] = new double*[nj];	//each data[i] points to an array of nj pointers-to-double
				for (int j=0;j<nj;j++)
				{	
					data[i][j] = new double[nk];	//each data[i][j] points to an array of nk doubles					
				}
			}
		 operator=(0);
		}
		
	//destructor, frees memory in reverse order
	~Field() {
		for (int i=0;i<ni;i++)
		{
			for (int j=0;j<nj;j++)
				delete data[i][j];	
			delete data[i];
		}
		
		delete[] data;	
	}
	
	double ** operator[] (int i) {return data[i];}

	 //overload the assignment operator
	 Field& operator= (double s) {
	 	for (int i=0;i<ni;i++)
	 		for (int j=0;j<nj;j++)
	 			for (int k=0;k<nk;k++)
	 				data[i][j][k] = s;
			return *this;
	 }

	const int ni,nj,nk;		//number of nodes
				
protected:
	double ***data;

};

#endif
