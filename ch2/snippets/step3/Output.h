#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "World.h"
#include "Species.h"

namespace Output {
	/*writes mesh data to a VTK image file*/
	void fields(World &world, std::vector<Species> &species);	
}

#endif
