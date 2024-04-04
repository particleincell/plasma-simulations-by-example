/*Field is a container for mesh node data
It contains functions for scatter/gather, and division by volume*/

#ifndef _FIELD_H
#define _FIELD_H

#include "World.h"
#include <fstream>
#include <mpi.h>

class Field 
{
public:
	
	/*constructor*/
	Field(World &world) :
		world(world)			/*sets this->world to world*/
	{
		data = new double**[world.ni];
		for (int i=0;i<world.ni;i++)
		{
			data[i] = new double*[world.nj];
			for (int j=0;j<world.nj;j++) data[i][j] = new double[world.nk];
		}		

		clear();
	}
	
	/*destructor, frees memory*/
	~Field()
	{
		for (int i=0;i<world.ni;i++)
		{
			for (int j=0;j<world.nj;j++) delete[] data[i][j];
			delete[] data[i];
		}
		delete[] data;
	}	
	
	/* scatters scalar value onto a field at logical coordinate lc*/
	void scatter(double lc[3], double value)
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		int k = (int)lc[2];
		double dk = lc[2]-k;
		
		data[i][j][k] += value*(1-di)*(1-dj)*(1-dk);
		data[i+1][j][k] += value*(di)*(1-dj)*(1-dk);
		data[i+1][j+1][k] += value*(di)*(dj)*(1-dk);
		data[i][j+1][k] += value*(1-di)*(dj)*(1-dk);
		data[i][j][k+1] += value*(1-di)*(1-dj)*(dk);
		data[i+1][j][k+1] += value*(di)*(1-dj)*(dk);
		data[i+1][j+1][k+1] += value*(di)*(dj)*(dk);
		data[i][j+1][k+1] += value*(1-di)*(dj)*(dk);		
	}

	/* gathers field value at logical coordinate lc*/
	double gather(double lc[3])
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		int k = (int)lc[2];
		double dk = lc[2]-k;
					
		/*gather electric field onto particle position*/
		double val = data[i][j][k]*(1-di)*(1-dj)*(1-dk)+
				data[i+1][j][k]*(di)*(1-dj)*(1-dk)+
				data[i+1][j+1][k]*(di)*(dj)*(1-dk)+
				data[i][j+1][k]*(1-di)*(dj)*(1-dk)+
				data[i][j][k+1]*(1-di)*(1-dj)*(dk)+
				data[i+1][j][k+1]*(di)*(1-dj)*(dk)+
				data[i+1][j+1][k+1]*(di)*(dj)*(dk)+
				data[i][j+1][k+1]*(1-di)*(dj)*(dk);
				
		return val;
	}
	
	/*sets all values to zero*/
	void clear()
	{
		for (int i=0;i<world.ni;i++)
			for (int j=0;j<world.nj;j++)
				for (int k=0;k<world.nk;k++) data[i][j][k] = 0;
	}
	
	/*shortcut to data*/
	double at(int i, int j, int k)
	{
		return data[i][j][k];
	}
	
	/*divides value by node volume*/
	void divideByVolume()
	{
		for (int i=0;i<world.ni;i++)
			for (int j=0;j<world.nj;j++)
				for (int k=0;k<world.nk;k++)
				{
					double dV=world.dh[0]*world.dh[1]*world.dh[2];
					/*divide be 2 at external boundaries*/
					if ((i==0 && world.mpi_i==0) ||
						(i==world.ni-1 && world.mpi_i==world.mpi_size_i-1)) dV*=0.5;
					if ((j==0 && world.mpi_j==0) ||
						(j==world.nj-1 && world.mpi_j==world.mpi_size_j-1)) dV*=0.5;
					if ((k==0 && world.mpi_k==0) ||
						(k==world.nk-1 && world.mpi_k==world.mpi_size_k-1)) dV*=0.5;

					data[i][j][k]/=dV;
				}			
	}

	/*updates boundaries across MPI processes*/
	void updateBoundaries()
	{
		const int TAG_DATA = 20;
		int ni=world.ni;	/*to save on typing*/
		int nj=world.nj;
		int nk=world.nk;

		for (int face=0;face<6;face++)
		{
			int target = world.getNeighborRank(face);
			/*update left*/
			if (target>=0)
			{
				int size;
				if (face/2==0) size=nj*nk;
				else if (face/2==1) size=ni*nk;
				else if (face/2==2) size=ni*nj;

				/*copy values along left boundary to send buffer*/
				double *send_buffer = new double[size];

				/*pack data, there is probably a more clever way of doing this*/
				double *s=send_buffer;	/*pointer to start*/
				switch (face)
				{
					case 0:  for (int j=0;j<nj;j++) for (int k=0;k<nk;k++) (*s++)=data[0][j][k];break;
					case 1:  for (int j=0;j<nj;j++) for (int k=0;k<nk;k++) (*s++)=data[ni-1][j][k];break;
					case 2:  for (int i=0;i<ni;i++) for (int k=0;k<nk;k++) (*s++)=data[i][0][k];break;
					case 3:  for (int i=0;i<ni;i++) for (int k=0;k<nk;k++) (*s++)=data[i][nj-1][k];break;
					case 4:  for (int i=0;i<ni;i++) for (int j=0;j<nj;j++) (*s++)=data[i][j][0];break;
					case 5:  for (int i=0;i<ni;i++) for (int j=0;j<nj;j++) (*s++)=data[i][j][nk-1];break;
				}

				double *recv_buffer = new double[size];

				/*transfer values*/
				MPI_Sendrecv(send_buffer, size, MPI_DOUBLE, target, TAG_DATA,
							 recv_buffer, size, MPI_DOUBLE, target, TAG_DATA,
							 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				/*add to our values*/
				double *r=recv_buffer;	/*pointer to start*/
				switch (face)
				{
					case 0:  for (int j=0;j<nj;j++) for (int k=0;k<nk;k++) data[0][j][k]+=(*r++);break;
					case 1:  for (int j=0;j<nj;j++) for (int k=0;k<nk;k++) data[ni-1][j][k]+=(*r++);break;
					case 2:  for (int i=0;i<ni;i++) for (int k=0;k<nk;k++) data[i][0][k]+=(*r++);break;
					case 3:  for (int i=0;i<ni;i++) for (int k=0;k<nk;k++) data[i][nj-1][k]+=(*r++);break;
					case 4:  for (int i=0;i<ni;i++) for (int j=0;j<nj;j++) data[i][j][0]+=(*r++);break;
					case 5:  for (int i=0;i<ni;i++) for (int j=0;j<nj;j++) data[i][j][nk-1]+=(*r++);break;
				}

				delete[] send_buffer;
				delete[] recv_buffer;
			}
		} /*for face*/

	}


	/*divides value by node volume*/
	void divideByField(Field *f)
	{
		for (int i=0;i<world.ni;i++)
			for (int j=0;j<world.nj;j++)
				for (int k=0;k<world.nk;k++)
				{
					double val = f->data[i][j][k];
					
					if (val!=0) data[i][j][k]/=val;
					else data[i][j][k]=0;
				}
			
	}

	/*writes out own data to file stream*/
	void write(std::ofstream &out)
	{	
		for (int k=0;k<world.nk;k++)
		{
			for (int j=0;j<world.nj;j++)
				for (int i=0;i<world.ni;i++) out<<data[i][j][k]<<" ";
			out<<"\n";
		}
	}

	int num_samples=0;  /*for averaged data*/
    
	/*add data to running average*/
    void updateAverage(Field *inst)
    {
        for (int i=0;i<world.ni;i++)
            for (int j=0;j<world.nj;j++)
                for (int k=0;k<world.nk;k++)
                {
                    data[i][j][k] = (num_samples*data[i][j][k]+inst->data[i][j][k])/(num_samples+1);
                }
        
        num_samples++;
    }
	
	double ***data;	/*data held by this field*/
	World &world;	/*this is a reference and needs to be set by the constructor*/
};


/*Helper class for 3D data*/
class Field3
{
public: 
	Field *f[3];
	
	Field3(World &world)
	{
		for (int i=0;i<3;i++) f[i]=new Field(world);
	}
	
	~Field3(){for (int i=0;i<3;i++) delete f[i];}

	/*sets all 3 components to zero*/
	void clear() {for (int i=0;i<3;i++) f[i]->clear();}
	
	/*scatter helper function*/
	void scatter(double lc[], double vec[], double s) 
	{
		for (int i=0;i<3;i++) f[i]->scatter(lc,vec[i]*s);
	}

	/*gather helper function*/
	void gather(double r[], double lc[]) 
	{
		for (int i=0;i<3;i++) r[i] = f[i]->gather(lc);
	}

	/*writes out data to a vtk file*/
	void write(std::ofstream &out)
	{
		World &world = f[0]->world;

		for (int k=0;k<world.nk;k++)
		{
			for (int j=0;j<world.nj;j++)
				for (int i=0;i<world.ni;i++) out<<f[0]->data[i][j][k]<<" "
												<<f[1]->data[i][j][k]<<" "
												<<f[2]->data[i][j][k]<<" ";
			out<<"\n";
		}
	}
};

#endif
