#ifndef _VOLUME_MESH_H
#define _VOLUME_MESH_H
#include <string>
#include <vector>
#include <map>
#include "Vec.h"

enum NodeType {NORMAL,OPEN,INLET,SPHERE};

/*definition of a node*/
struct Node
{
	Node(double x, double y, double z): pos{x,y,z} {}
	double3 pos;	/*node position*/
	int type = NORMAL;
};

/*definition of a triangular surface element*/
class VolumeMesh;		//forward definition
struct Tri
{
	Tri (int n1, int n2, int n3): con{n1,n2,n3} {}
	int con[3];
	double area = 0;
	double3 normal; //normal vector
	double3 randomPos(const VolumeMesh &vm); 	//returns a random point on the triangle
};

/*definition of a tetrahedral volume element*/
struct Tet
{
	Tet (int n1, int n2, int n3, int n4): con{n1,n2,n3,n4} {}
	int con[4];
	double volume = 0;	//cell volume

	/*data structures to hold precomputed 3x3 determinants*/
	double alpha[4], beta[4], gamma[4], delta[4];

	/*cell neighbors*/
	int neighbor[4]; /*cell id attached to the face opposite the i-th node*/
};


struct Group{
	std::vector<int> elements;
};

/*defines the logical coordinate*/
struct LCord {
	int cell_id = 0;				//cell id
	double L[4] = {0,0,0,0};	//shape factors
	LCord (int c=0):cell_id{c} {}
	explicit operator bool() {return cell_id>=0;} //if(lc) checks for cell_id>=0
};

/*definition of a volume*/
class VolumeMesh
{
public:
	std::vector <Node> nodes;
	std::vector <Tri> tris;
	std::vector <Tet> tets;
	std::map<std::string,Group> groups;

	VolumeMesh(std::string file_name) {loadUNV(file_name);}

	/*computes logical coordinate*/
	LCord XtoL(const double3 &pos, int e, bool search=true);

protected:
	bool loadUNV(std::string file_name);
	void loadNodes(std::ifstream &in);
	void loadElements(std::ifstream &in);
	void loadGroups(std::ifstream &in);
	void skipDataset(std::ifstream &in);
	void init();
	std::map<int,int> lut_nodes; //look up table
	std::map<int,int> lut_tris;
	std::map<int,int> lut_tets;
};

#endif
