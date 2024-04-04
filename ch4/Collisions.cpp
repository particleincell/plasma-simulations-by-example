#include <vector>
#include <math.h>
#include <iostream>
#include "Collisions.h"

using namespace std;


//approximates ionization, not mass conserving
void ChemistryIonize::apply(double dt) {
	double dV = world.getCellVolume();
	//loop over all cells
	for (int i=0;i<world.ni-1;i++)
		for (int j=0;j<world.nj-1;j++)
			for (int k=0;k<world.nk-1;k++) {
				//evaluate electron and neutral density at cell center
				double na = neutrals.den.gather({i+0.5,j+0.5,k+0.5});
				double ne = world.rho.gather({i+0.5,j+0.5,k+0.5})/Const::QE; //assume QN
				double dni = rate*na*ne*dt;

				/*number of macroparticles to create*/
				int num_p = (int)(dni*dV/ions.mpw0 + rnd());
				for (int p=0;p<num_p;p++) {
					/*sample a random particle in the cell*/
					double3 pos = world.pos(i+rnd(),j+rnd(),k+rnd());
					double3 vel {0,0,0}; //should sample from Maxwellian instead
					ions.addParticle(pos,vel,ions.mpw0);
				}

			}
}

//return 1e-18*pow(rel_g,-0.5);

/*contant cross-section for the demo*/
double evalSigma(double rel_g)
{
	return 1e-16;
}


/*MCC*/
void MCC_CEX::apply(double dt)
{
		/*set pointers to target data*/
		Field &target_den = target.den;
		Field &target_temp = target.T;
		Field3 &target_vel = target.vel;

		/*loop over all particles*/
        for (Particle &part:source.particles)
        {
			/*get target velocity*/
			double3 lc = world.XtoL(part.pos);
            			
			double3 vt = target_vel.gather(lc);

            /*get target density*/
            double nn = target_den.gather(lc);
           
			/*compute cross-section, function of relative velocity*/
            double3 v_rel = part.vel-vt;
			double v_rel_mag = mag(v_rel);
            
            /*evaluate cross-section */
            double sigma = 1e-16;
            
            /*compute probability*/
            double P = 1 - exp(-nn*sigma*v_rel_mag*dt);
            
            /*compare to a random number to see if collision happened*/
            if (P>=rnd())
            {
                /*sample a virtual target particle*/
                double T_target = target_temp.gather(lc);	//sample target temperature
                double3 v_target = target.sampleIsotropicVel(T_target);
                //part.vel = v_target; 	//CEX collision
                part.vel = 0;
            }
        }
    }
    

/* collides two particles*/
void DSMC_MEX::collide(double3 &vel1, double3 &vel2, double mass1, double mass2)
{
	double3 cm = (mass1*vel1 + mass2*vel2)/(mass1+mass2);

	/*relative velocity, magnitude remains constant through the collision*/
	double3 cr = vel1 - vel2;
	double cr_mag = mag(cr);

	/*pick two random angles, per Bird's VHS method*/
	double cos_chi = 2*rnd()-1;
	double sin_chi = sqrt(1-cos_chi*cos_chi);
	double eps = 2*Const::PI*rnd();

	/*perform rotation*/
	cr[0] = cr_mag*cos_chi;
	cr[1] = cr_mag*sin_chi*cos(eps);
	cr[2] = cr_mag*sin_chi*sin(eps);

	/*post collision velocities*/
	vel1 = cm + mass2/(mass1+mass2)*cr;
	vel2 = cm - mass1/(mass1+mass2)*cr;
}


/*DSMC*/
void DSMC_MEX::apply(double dt)
{
	/*first we need to sort particles to cells*/
	vector<Particle*> *parts_in_cell;
	int n_cells = (world.ni-1)*(world.nj-1)*(world.nk-1);
	parts_in_cell = new vector<Particle*> [n_cells];

	/*sort particles to cells*/
	for (Particle &part:species.particles)
	{
		int c = world.XtoC(part.pos);
		parts_in_cell[c].push_back(&part);
	}

	double sigma_cr_max_temp = 0;	/*reset for max computation*/
	double dV = world.getCellVolume();	/*internal cell volume*/
	double Fn = species.mpw0;	/*specific weight, using Bird's notation*/
	int num_cols=0;	/*reset collision counter*/
					
	/*now perform collisions*/
	for (int c=0;c<n_cells;c++)
	{
		vector<Particle*> &parts = parts_in_cell[c];
		int np = parts.size();
		if (np<2) continue;

		/*compute number of groups according to NTC*/
		double ng_f = 0.5*np*np*Fn*sigma_cr_max*dt/dV;
		int ng = (int)(ng_f+0.5);	/*number of groups, round*/
	
		/*assumes at least two particles per cell*/
		for (int g=0;g<ng;g++)
		{
			int p1, p2;
			p1 = (int)(rnd()*np);		/*returns some number between 0 and np-1 inclusive*/
		
			do {
				p2 = (int)(rnd()*np);
			} while (p2==p1);

			/*compute relative velocity*/
			double3 cr_vec = parts[p1]->vel - parts[p2]->vel;
			double cr = mag(cr_vec);

			/*evaluate cross section*/
			double sigma = evalSigma(cr);

			/*eval sigma_cr*/
			double sigma_cr=sigma*cr;

			/*update sigma_cr_max*/
			if (sigma_cr>sigma_cr_max_temp)
				sigma_cr_max_temp=sigma_cr;

			/*eval prob*/
			double P=sigma_cr/sigma_cr_max;

			/*did the collision occur?*/
			if (P>rnd())
			{
				num_cols++;
				collide(parts[p1]->vel,parts[p2]->vel,species.mass, species.mass);
			}
		}
	}

	delete[] parts_in_cell;

	if (num_cols){
		sigma_cr_max = sigma_cr_max_temp;
	}

}
