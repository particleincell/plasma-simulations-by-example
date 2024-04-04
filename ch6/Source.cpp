#include <iostream>
#include <string>
#include <stdexcept>
#include "Source.h"
#include "World.h"

using namespace std;

ColdBeamSource::ColdBeamSource(Species &species, World &world, double v_drift, double den, Group &group) :
sp{species}, world{world}, v_drift{v_drift}, den{den}, src_group{group} {
	//compute total source area
	A_tot = 0;
	for (int e:group.elements)
		A_tot += world.vm.tris[e].area;
}


//samples monoenergetic particles according to a prescribed density
void ColdBeamSource::sample()
{
	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	double num_real = den*v_drift*A_tot*world.getDt();

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	VolumeMesh &vm = world.vm;

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		Tri *tri;
		do {
			//pick random entry in the group
			int ge = (int)(rnd()*src_group.elements.size());
			tri = &vm.tris[src_group.elements[ge]];
		} while (tri->area/A_tot<rnd());

		double3 pos = tri->randomPos(vm);
		double3 vel = tri->normal*v_drift;
		sp.addParticle(pos,vel,sp.mpw0);
	}

}


