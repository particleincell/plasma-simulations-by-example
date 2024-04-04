#ifndef _WORLD_H
#define _WORLD_H

#include "Field.h"

//World.h
class World {
public:
  World (int ni, int nj, int nk);  //constructor
  
  //sets the mesh span, also recomputes cell spacing
  void setExtents(const double3 x0, const double3 xm);

  double3 getX0() const {return double3(x0);}
  double3 getXm() const {return double3(xm);}
  double3 getXc() const {return double3(xc);}
  double3 getDh() const {return double3(dh);}
    
  const int nn[3];	 	//number of nodes  
  const int ni,nj,nk;	//number of nodes in individual variables
  
  Field phi;		//potential
  Field rho;		//charge density
  Field node_vol;	//node volumes
  Field3 ef;		//electric field

protected:
  double x0[3];	 	//mesh origin
  double dh[3];	 	//cell spacing
  double xm[3];		//mesh max bound
  double xc[3];		//domain centroid
};

#endif
