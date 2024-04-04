//Main.cpp
#include "World.h"
#include "Field.h"

int main(/*...*/)
{
	World world(21,21,21);		//calls World constructor to create a variable 'world'
	world.setExtents(-0.1,-0.1,-0.1,0.1,0.1,0.2); //x1=(-0.1,-0.1,0), x2=(0.1,0.1,0.2)

	Field phi(21,21,21);		//calls Field constructor
	phi = 0;
	phi[10][2][3]=3;
}

