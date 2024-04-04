#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>



using namespace std;

int main() {

	double theta1 = 0;
	double theta2 = 270;
	int nk = 1;		//theta
	int ni = 6;
	int nj = 4;
	double L = 1;
	double r = 0.4;
	const double PI = acos(-1.0);

	ofstream out("plane.vts");

	out<<"<VTKFile type=\"StructuredGrid\" byte_order=\"LittleEndian\">"<<endl;
	out<<"<StructuredGrid WholeExtent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 "<<nk-1<<"\">"<<endl;
	out<<"<Piece Extent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 "<<nk-1<<"\">"<<endl;
	out<<"<Points>"<<endl;
	out<<"<DataArray Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">"<<endl;

	for (int k=0;k<nk;k++) 
	{
		double theta = (nk>1)?(PI/180.0)*(theta1+k/(nk-1.0)*(theta2-theta1)):0;
		for (int j=0;j<nj;j++)		
			for (int i=0;i<ni;i++) 
			{
				double x = i/(ni-1.0)*L;
				double y = j/(nj-1.0)*r*cos(theta);
				double z = j/(nj-1.0)*r*sin(theta);					
				out<<x<<" "<<y<<" "<<z<<endl;
			}
	}

	out<<"</DataArray>"<<endl;
	out<<"</Points>"<<endl;
	out<<"</Piece>"<<endl;
	out<<"</StructuredGrid>"<<endl;
	out<<"</VTKFile>"<<endl;
	return 0;
}

