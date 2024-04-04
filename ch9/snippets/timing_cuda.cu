#include <iostream>
#include <chrono>
#include <cuda_runtime.h>
using namespace std;

struct Particle {
	float pos[3] = {0,0,0};
	float vel[3] = {0,0,0};
};

//kernel code to run on the GPU
__global__ void gpu_push(Particle *particles, float dt, size_t N) {
   int p = blockIdx.x*blockDim.x + threadIdx.x;
   if (p<N) {
	   Particle &part = particles[p];
	   for (int i=0;i<3;i++)
		part.pos[i] += part.vel[i]*dt;
	}
}

//code to push a single particle
void push(Particle *part, double dt) {
	for (int i=0;i<3;i++)
		part->pos[i] += part->vel[i]*dt;
}

int main(int n_args, char *args[]) {

  cudaFree(0);
  size_t num_particles = 1000000;
 //Particle *particles = new Particle[num_particles];
  Particle *particles;
  cudaHostAlloc(&particles,sizeof(Particle)*num_particles,cudaHostAllocDefault); //allocate pinned memory

  //set some initial values
  for (size_t i=0;i<num_particles;i++)
	particles[i].vel[0]=1/(double)num_particles;

  const float dt = 0.1;

  //*** CPU particle push ***
   auto start_cpu = chrono::system_clock::now();
   for (size_t i=0;i<num_particles;i++) push(&particles[i],dt);
   auto end_cpu = chrono::system_clock::now();

   //*** GPU particle push ***
  auto start_gpu = chrono::system_clock::now();
  Particle *devParticles;
  cudaMalloc((void**)&devParticles, sizeof(Particle)*num_particles);
  cudaMemcpy(devParticles,particles,sizeof(Particle)*num_particles,cudaMemcpyHostToDevice);

  const int threads_per_block = 1024;
  int num_blocks = (num_particles-1)/threads_per_block + 1;
  cout<<"Creating "<<num_blocks*threads_per_block<<" threads"<<endl;
  gpu_push<<<num_blocks,threads_per_block>>>(devParticles, dt, num_particles);
  cudaMemcpy(particles,devParticles,sizeof(Particle)*num_particles,cudaMemcpyDeviceToHost);

  auto end_gpu = chrono::system_clock::now();
    
  //output timing info
  std::chrono::duration<double,std::nano> elapsed_cpu = end_cpu - start_cpu;
  std::chrono::duration<double,std::nano> elapsed_gpu = end_gpu - start_gpu;
  cout<<"Time per particle on CPU: "<<elapsed_cpu.count()/num_particles<<" (ns)"<<endl;
  cout<<"Time per particle on GPU: "<<elapsed_gpu.count()/num_particles<<" (ns)"<<endl;

//  delete[] particles;
  cudaFreeHost(particles);
 
  return 0;
}
