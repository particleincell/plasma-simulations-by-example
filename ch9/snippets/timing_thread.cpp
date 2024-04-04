#include <thread>
#include <iostream>
#include <vector>
#include <chrono>
#include <sstream>

using namespace std;
using dvector = vector<double>;

struct Particle {
	double pos[3] = {0,0,0};
	double vel[3] = {0,0,0};
};

//code to push a single particle
void push(Particle *part, double dt) {
	for (int i=0;i<3;i++)
		part->pos[i] += part->vel[i]*dt;
}

int main(int n_args, char *args[]) {

  size_t num_particles = 100000;
  vector<Particle> particles(num_particles);
  //set some initial values
  for (size_t i=0;i<num_particles;i++)
	particles[i].vel[0]=1/(double)num_particles;

  const double dt = 0.1;
  //*** particle push timing ***
  auto start_push = chrono::system_clock::now();  
  for (size_t i=0;i<num_particles;i++)
	push(&particles[i],dt);
  auto end_push = chrono::system_clock::now();  

  //*** thread timing ***
  auto start_threads = chrono::system_clock::now();  
  for (int i=0;i<num_particles;i++)
  {
    thread t(push,&particles[i],dt);
	t.join(); 
  }
  auto end_threads = chrono::system_clock::now();
    
  //output timing info
  std::chrono::duration<double,std::micro> elapsed_push = end_push - start_push;
  std::chrono::duration<double,std::micro> elapsed_threads = end_threads - start_threads;
  cout<<"Time to push one particle: "<<elapsed_push.count()/num_particles<<" (us)"<<endl;
  cout<<"Time to create a thread: "<<(elapsed_threads.count()-elapsed_push.count())/num_particles<<" (us)"<<endl;
 
  return 0;
}
