//ch1_01.cpp
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;
using dvector = vector<double>;

/*function prototypes*/
bool outputCSV(double x0, double dx, const dvector &phi, const dvector &rho, const dvector &ef);

/*constants*/
namespace Const
{
	const double QE = 1.602176565e-19;	// C, electron charge
	const double EPS_0 = 8.85418782e-12;// C/V/m, vacuum permittivity
	const double ME = 9.10938215e-31;	// kg, electron mass
};

using namespace Const;	//to avoid having to write Const::QE

/*main*/
int main()
{
	const int ni = 21;	//number of nodes
	const double x0 = 0;	//origin
	const double xd = 0.1;	//opposite end	
	double dx = (xd-x0)/(ni-1);	//node spacing

	dvector phi(ni);
	dvector rho(ni,QE*1e12);
	dvector ef(ni);
	
	//ouput to a CSV file for plotting
	outputCSV(x0,dx,phi,rho,ef);	
	
	return 0;	//normal exit
}

/*outputs the given fields to a CSV file, returns true if ok*/
bool outputCSV(double x0, double dx, const dvector &phi, const dvector &rho, const dvector &ef)
{
	ofstream out("results.csv");	//open file for writing
	if (!out)
	{
		cerr<<"Could not open output file!"<<endl; 
		return false;
	}
	
	out<<"x,phi,rho,ef\n";		//write header
	for (int i=0;i<phi.size();i++)
	{
		out<<x0+i*dx; //write i-th position
		out<<","<<phi[i]<<","<< rho[i]<<","<<ef[i]; //write values
		out<<"\n";	//new line, not using endl to avoid buffer flush
	}
	
	//file closed automatically when "out" variable is destroyed
	return true;
}
