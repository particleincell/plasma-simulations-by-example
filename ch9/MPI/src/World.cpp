/*defines the simulation domain*/
#include "World.h"
#include "Field.h"
#include <math.h>
	
/*constructor*/
World::World(int ni, int nj, int nk)	
{
	this->ni=ni;this->nj=nj;this->nk=nk;

	/*allocate memory*/
	object = new int**[ni];
	for (int i=0;i<ni;i++)
	{
		object[i] = new int*[nj];
		for (int j=0;j<nj;j++) 
		{
			object[i][j] = new int[nk];
			for (int k=0;k<nk;k++) object[i][j][k] = 0;
		}
	}		


	phi = new Field(*this);
	rhoi = new Field(*this);
	ef = new Field3(*this);
}

/*destructor, frees memory*/
World::~World()
{
	for (int i=0;i<ni;i++)
	{
		for (int j=0;j<nj;j++) delete[] object[i][j];
		delete[] object[i];
	}
	delete[] object;

	delete phi;
	delete rhoi;
	delete ef;
}

/*marks k=0 plane as dirichlet*/
void World::AddInlet()
{
	/*only on first z domain*/
	if (mpi_k!=0) return;

    for (int i=0;i<ni;i++)
        for (int j=0;j<nj;j++)
        {
            object[i][j][0] = 1;
            phi->data[i][j][0] = 0;	/*0V*/
        }
}
	
/*sugarcubes a sphere centered at (x0,y0,z0)*/
void World::AddSphere(double x0, double y0, double z0, double radius, double phi_sphere)
{
    /*save sphere parameters*/
    sphere_x0[0] = x0;
    sphere_x0[1] = y0;
    sphere_x0[2] = z0;
    sphere_rad = radius;

    for (int i=0;i<ni;i++)
        for (int j=0;j<nj;j++)
            for (int k=0;k<nk;k++)
            {
                /*compute node position*/
                double r[3];
				pos(r,i,j,k);
                if (inSphere(r))
                {
                    object[i][j][k] = 2;
                    phi->data[i][j][k] = phi_sphere;
                }
            }
}
	
/*returns true if point r is inside or on the sphere*/
bool World::inSphere(double r[3])
{
    double dx = r[0]-sphere_x0[0];
    double dy = r[1]-sphere_x0[1];
    double dz = r[2]-sphere_x0[2];

    double r_mag = sqrt(dx*dx + dy*dy + dz*dz);
    if (r_mag<=sphere_rad) return true;
    return false;
}

