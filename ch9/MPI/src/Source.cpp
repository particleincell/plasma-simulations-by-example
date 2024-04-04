#include "Source.h"

double num_rem;		/*remainder number of ions*/	
/*samples ions at the k=0 edge such that local plasma density equals PLASMA_DEN*/
void Source::sample()
{
	/*only sample on global k=0 face*/
	if (world.mpi_k!=0) return;

    /*compute area of the k=0 face*/
    double area = (world.ni-1)*world.dh[0]*(world.nj-1)*world.dh[1];

    /*number of real ions per sec, given prescribed density and velocity*/
    double num_per_sec = den0*v_drift*area;

    /*number of ions to generate in this time step*/
    double num_real = num_per_sec*world.dt;

    /*fraction number of macroparticles*/
    double fnum_mp = num_real/species.spwt + num_rem;

    /*integer number of macroparticles*/
    int num_mp = (int)fnum_mp;

    /*update remainder*/
    num_rem = fnum_mp - num_mp;

    /*sample particles*/
    for (int p=0;p<num_mp;p++)
    {
        /*sample random position on the inlet face*/
        double pos[3];
        pos[0] = world.x0[0]+rnd()*(world.ni-1)*world.dh[0];
        pos[1] = world.x0[1]+rnd()*(world.nj-1)*world.dh[1];
        pos[2] = world.x0[2];

        /*load isotropic thermal distribution, 
        http://www.particleincell.com/blog/2012/isotropic-velocity/*/
        /*pick a random angle*/
		double theta = 2*PI*rnd();
 
        double R = -1.0+2*rnd();    /*pick a random direction for n[2]*/
        double a = sqrt(1-R*R);
 
	    double n[3];
        n[0] = cos(theta)*a;
        n[1] = sin(theta)*a;
        n[2] = R;
            
        double vel[3];
        vel[0] = v_th*n[0];
        vel[1] = v_th*n[1];
        vel[2] = v_th*n[2] + v_drift;
			
        /*reverse if going in wrong direction*/
        if (vel[2]<0) vel[2]=-vel[2];

        /*rewind velocity*/
        double lc[3];
		world.XtoL(lc, pos);

        double ef_part[3];

        world.ef->gather(ef_part,lc);

        for (int i=0;i<3;i++)
            vel[i] -= species.charge/species.mass*ef_part[i]*(0.5*world.dt);

        /*add to list*/
        species.addParticle(pos,vel);						
    }		
}
