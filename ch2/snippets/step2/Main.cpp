//Main.cpp
#include "World.h"
#include "Output.h"

int main() {
  /*initialize domain*/
  World world(21,21,21);
  world.setExtents({-0.1,-0.1,0.0},{0.1,0.1,0.2});
  
  Output::fields(world);
  return 0;
}


