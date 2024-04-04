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
	MCC_CEX(KineticSpecies &source, Species &target, World &world) :
		source{source}, target{target}, world{world} {}
	void apply(double dt);
protected:
	KineticSpecies &source;
	Species &target;
	World &world;
};


#endif
