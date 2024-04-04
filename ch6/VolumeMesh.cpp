#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include "VolumeMesh.h"
#include "FESolver.h"

using namespace std;

/*loads volume mesh from a .unv file*/
bool VolumeMesh::loadUNV(string file_name) {
	ifstream in(file_name);
	if (!in) {cerr<<"Failed to open "<<file_name<<endl;return false;}

	while (in) {
		string line;
		getline(in,line);
		if (in.eof()) break;
		if (line != "    -1") {cerr<<"Data file mismatch"<<endl;return false;}
		getline(in,line);
		stringstream ss(line);
		int code;	ss>>code;

		switch (code) {
			case 2411: loadNodes(in); break;
			case 2412: loadElements(in); break;
			case 2467: loadGroups(in); break;
			default: skipDataset(in); break;
		}
	}

	init();  //compute node volumes
	return true;
}

/*loads nodes from dataset 2411*/
void VolumeMesh::loadNodes(ifstream &in) {
	string line;
	const double mm = 1e-3;	 //data is in mm

	while (in) {
		getline(in, line);
		if (line == "    -1") break;
		stringstream ss(line);		//attach stream to the string
		int id;	ss>>id;	//read node_id
		getline(in,line);	//read next line
		stringstream ss2(line);
		double x,y,z;
		ss2>>x>>y>>z;
		lut_nodes[id] = nodes.size();	//add to look up table
		nodes.emplace_back(x*mm,y*mm,z*mm);
	}
	cout<<"Loaded "<<nodes.size()<<" nodes"<<endl;
}

/*loads elements from dataset 2412*/
void VolumeMesh::loadElements(ifstream &in) {
	string line;
	while (in) {
		getline(in, line);
		if (line  == "    -1") break;
		stringstream ss(line);
		int id, temp, num_nodes;
		ss>>id>>temp>>temp>>temp>>temp>>num_nodes;

		if (num_nodes==2) getline(in, line);	//skip beam info
		getline(in, line);
		stringstream ss2(line);
		int n1,n2,n3,n4;
		//load  nodes
		if (num_nodes==3) {
			ss2>>n1>>n2>>n3;
			lut_tris[id] = tris.size();
			tris.emplace_back(lut_nodes[n1],lut_nodes[n2],lut_nodes[n3]);
		}
		else if (num_nodes==4) {
			ss2>>n1>>n2>>n3>>n4;
			lut_tets[id] = tets.size();
			tets.emplace_back(lut_nodes[n1],lut_nodes[n2],lut_nodes[n3],lut_nodes[n4]);
		}
		else {
			/*skip non 3 and 4 node elements*/
		}
	}
	cout<<"Loaded "<<tris.size()<<" triangles and "<<tets.size()<<" tetrahedrons"<<endl;
}

/*loads groups from dataset 2467*/
void VolumeMesh::loadGroups(ifstream &in) {
	string line;
	while (in) {
		getline(in, line);
		if (line  == "    -1")
			break;

		//get number of elements in this group
		int temp,num_eles;
		stringstream ss(line);
		ss>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>num_eles;
		string name;
		getline(in,name);		//get group name
		Group group;
		for (int i=0;i<(num_eles+1)/2;i++) /*number of lines*/
		{
			getline(in,line);
			stringstream ss2(line);
			int type, id, temp;
			ss2>>type>>id>>temp>>temp;
			if (type==8) {
				//assume surface element groups
				group.elements.emplace_back(lut_tris[id]);
			}
			ss2>>type>>id>>temp>>temp;	//second item
			if (type==8) {
				//assume surface element groups
				group.elements.emplace_back(lut_tris[id]);
			}
		}
		//add to list
		groups[name] = group;
		cout<<"Loaded group "<<name<<endl;
	}
}

/*skips data from an arbitrary dataset until the terminating -1*/
void VolumeMesh::skipDataset(ifstream &in) {
	string line;
	while (true && in) {
		getline(in, line);
		if (line=="    -1") break;
	}
}

/*loads and initializes volume mesh*/
void VolumeMesh::init()
{
	/*fix Dirichlet nodes by looping over map*/
	for (auto &m:groups) {
		NodeType type;
		if (m.first=="inlet") type = NodeType::INLET;
		else if (m.first=="sphere") type = NodeType::SPHERE;
		else {cerr<<"Ignoring unknown group "<<m.first<<endl;continue;}

		Group &group = m.second;
		for (int &i:group.elements) //loop over all triangles in the group
			for (int j=0;j<3;j++) //loop over nodes
				nodes[tris[i].con[j]].type = type;
	}

	/*compute cell volumes*/
	for (Tet &tet: tets)
	{
		double M[4][4];

		/*set first column to 1*/
		for (int i=0;i<4;i++) M[i][0] = 1;

		/*loop over vertices*/
		for (int v=0;v<4;v++)
		{
			for (int dim=0;dim<3;dim++)
			{
				M[0][dim+1] = nodes[tet.con[0]].pos[dim];
				M[1][dim+1] = nodes[tet.con[1]].pos[dim];
				M[2][dim+1] = nodes[tet.con[2]].pos[dim];
				M[3][dim+1] = nodes[tet.con[3]].pos[dim];
			}
		}

		/*volume is (1/6)*det4(M)*/
		tet.volume = (1.0/6.0)*utils::det4(M);

		/*flip 0123 to 0321 if negative volume*/
		if (tet.volume<0) {int t=tet.con[1];tet.con[1]=tet.con[3];tet.con[3]=t;tet.volume=-tet.volume;}
	}

	/*precompute 3x3 determinants for LC computation, needs to be done after connectivity flip*/
	for (Tet &tet:tets)
	{
		double M[3][3];
		/*loop over vertices*/
		for (int v=0;v<4;v++)
		{
			int v1,v2,v3;
			switch (v)	{
				case 0: v1=1;v2=2;v3=3;break; //Vp123
				case 1: v1=3;v2=2;v3=0;break; //Vp320
				case 2: v1=3;v2=0;v3=1;break; //Vp301
				case 3: v1=1;v2=0;v3=2;break; //Vp102
			}

			double3 p1 = nodes[tet.con[v1]].pos;
			double3 p2 = nodes[tet.con[v2]].pos;
			double3 p3 = nodes[tet.con[v3]].pos;

			/*alpha*/
			M[0][0] = p1[0]; M[0][1] = p1[1]; M[0][2] = p1[2];
			M[1][0] = p2[0]; M[1][1] = p2[1]; M[1][2] = p2[2];
			M[2][0] = p3[0]; M[2][1] = p3[1]; M[2][2] = p3[2];
			tet.alpha[v] = utils::det3(M);

			/*beta*/
			M[0][0] = 1; M[0][1] = p1[1]; M[0][2] = p1[2];
			M[1][0] = 1; M[1][1] = p2[1]; M[1][2] = p2[2];
			M[2][0] = 1; M[2][1] = p3[1]; M[2][2] = p3[2];
			tet.beta[v] = utils::det3(M);

			/*gamma*/
			M[0][0] = 1; M[0][1] = p1[0]; M[0][2] = p1[2];
			M[1][0] = 1; M[1][1] = p2[0]; M[1][2] = p2[2];
			M[2][0] = 1; M[2][1] = p3[0]; M[2][2] = p3[2];
			tet.gamma[v] = utils::det3(M);

			/*delta*/
			M[0][0] = 1; M[0][1] = p1[0]; M[0][2] = p1[1];
			M[1][0] = 1; M[1][1] = p2[0]; M[1][2] = p2[1];
			M[2][0] = 1; M[2][1] = p3[0]; M[2][2] = p3[1];
			tet.delta[v] = utils::det3(M);
		}
	}


	/*set surface areas and normals*/
	for (Tri &tri:tris) {
		double3 v0 = nodes[tri.con[0]].pos;
		double3 v1 = nodes[tri.con[1]].pos;
		double3 v2 = nodes[tri.con[2]].pos;
		double3 c = cross(v1-v0,v2-v0);
		tri.area = 0.5*mag(c);
		tri.normal = unit(c);	//unit vector
	}

	/*build cell connectivity, there is probably a faster way*/
	cout<<"Building cell connectivity"<<endl;

	/*reset cell neighbors*/
	for (Tet &tet:tets) { for (int i=0;i<4;i++) tet.neighbor[i]=-1;	}

	/*set cell neighbors*/
	for (size_t l=0;l<tets.size();l++)
	{
		Tet &tet = tets[l];
		int v1,v2,v3;
		for (int v=0;v<4;v++)
		{
			/*skip if already set*/
			if (tet.neighbor[v]>=0) continue;

			switch(v)
			{
				case 0: v1=1;v2=2;v3=3;break;
				case 1: v1=2;v2=3;v3=0;break;
				case 2: v1=3;v2=0;v3=1;break;
				case 3: v1=0;v2=1;v3=2;break;
			}

			/*loop over the remaining tets looking for one with these three vertices*/
			for (size_t m=l+1;m<tets.size();m++)
			{
				Tet &other = tets[m];

				bool matches[4] = {false,false,false,false};
				int count = 0;
				for (int k=0;k<4;k++)
				{
					if (other.con[k]==tet.con[v1] ||
					    other.con[k]==tet.con[v2] ||
						other.con[k]==tet.con[v3]) {count++;matches[k]=true;}
				}

				/*if three vertices match*/
				if (count==3)
				{
					tet.neighbor[v] = m;

					/*set the cell connectivity for the index without a matching vertex to l*/
					for (int k=0;k<4;k++)
						if(!matches[k]) other.neighbor[k] = l;
				}
			}
		}
	}

	/*mark nodes on open faces as open*/
	for (size_t e=0;e<tets.size();e++)
	{
		Tet &tet = tets[e];
		for (int v=0;v<4;v++)
			if (tet.neighbor[v]<0)	/*no neighbor*/
			{
				for (int i=0;i<4;i++)
				{
					if (i!=v) nodes[tet.con[i]].type=OPEN;
				}
			}
	}
}

/*converts physical coordinate to logical, returns LCord(-1) if search failed*/
LCord VolumeMesh::XtoL(const double3 &pos, int e, bool search)
{
	/*first try the current tetrahedron*/
	Tet &tet = tets[e];
	LCord lc(e);

	bool found = true;
	/*loop over vertices*/
	for (int i=0;i<4;i++)
	{
		lc.L[i] = (1.0/6.0)*(tet.alpha[i] - pos(0)*tet.beta[i] +
				pos(1)*tet.gamma[i] - pos(2)*tet.delta[i])/tet.volume;
		if (lc.L[i]<0 || lc.L[i]>1.0) found=false;
	}

	if (found) return lc;
	if (!search) return LCord(-1);

	/*we are outside the last known tet, find most negative weight*/
	int min_i = 0;
	double min_lc = lc.L[0];
	for (int i=1;i<4;i++)
		if (lc.L[i]<min_lc) {min_lc=lc.L[i];min_i=i;}

	/*is there a neighbor in this direction?*/
	if (tet.neighbor[min_i]>=0)
		return XtoL(pos,tet.neighbor[min_i]);

	/*no neighbor in this direction*/
	return LCord(-1);
}



/** returns random position in an element,
    quads are dealt with by picking a random sub-triangle
	based on http://www.cs.princeton.edu/~funk/tog02.pdf
    and https://stackoverflow.com/questions/4778147/sample-random-point-in-triangle
  */
double3 Tri::randomPos(const VolumeMesh &vm)
{
	//triangle vertex positions
	double3 verts[3] = {vm.nodes[con[0]].pos, vm.nodes[con[1]].pos, vm.nodes[con[2]].pos};

	double eps = 0.0001;	//avoid sampling on edges, 0 implies all the way to the edge
	double r1 = eps + rnd()*(1-2*eps);
    double r2 = eps + rnd()*(1-2*eps);
	double fac0 = (1-sqrt(r1));
	double fac1 = sqrt(r1)*(1-r2);
	double fac2 = sqrt(r1)*r2;
	return fac0*verts[0] + fac1*verts[1] + fac2*verts[2];
}

