#ifndef _COLLISIONS_H
#define _COLLISIONS_H

#include <vector>
#include "World.h"
#include "Species.h"
#include "Field.h"

class Sigma {
public:
	virtual double operator() (double g) {return 0;};
	virtual ~Sigma () {};
};

class SigmaPoly: public Sigma {
	SigmaPoly (const std::vector<double> &coeffs) : coeffs{coeffs} {}

	double operator() (double g) {
		double r = coeffs[0];
		for (size_t i=1;i<coeffs.size();i++) {
			r+=coeffs[i]*g; g*=g;
		}
		return r;
	}


protected:
		std::vector<double> coeffs;
};



class Interaction {
public:
	virtual void apply(double dt) = 0;
	virtual ~Interaction() {}
};

//MCC Charge Exchange Collision
class MCC_CEX: public Interaction {
public:
	MCC_CEX(Species &source, Species &target, World &world) :
		source{source}, target{target}, world{world} {}
	void apply(double dt);
protected:
	Species &source;
	Species &target;
	World &world;
};

class DSMC_MEX: public Interaction {
public:
	DSMC_MEX(Species &species, World &world) : species{species}, world{world} {

		mr = species.mass*species.mass/(species.mass + species.mass);
		c[0] = 4.07e-10;
		c[1] = 0.77;
		c[2]= 2*Const::K*273.15/mr;	//Bird's reference params at 273.15 K
		c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
	}
	void apply(double dt);
protected:
	void collide(double3 &vel1, double3 &vel2, double mass1, double mass2);
	double evalSigma(double g_rel) {
		return Const::PI*c[0]*c[0]*pow(c[2]/(g_rel*g_rel),c[1]-0.5)/c[3];

	}

	double sigma_cr_max = 1e-14;	//some initial value
	double mr;
	double c[4];
	Species &species;
	World &world;
};

class ChemistryIonize: public Interaction {
public:
	ChemistryIonize(Species &neutrals, Species &ions, World &world, double rate):
	neutrals{neutrals},ions{ions},world{world},rate{rate} { }
	void apply(double dt);
protected:
	Species &neutrals;
	Species &ions;
	World &world;
	double rate;
};


#endif
