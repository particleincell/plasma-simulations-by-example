#include <cuda_runtime.h>
#include <iostream>
using namespace std;

//code to run on the GPU
__global__ void add(float a, float b, float *c) {
   *c = a + b;
}

int main(int n_args, char *args[]) 
{
	float *dev_c;
	cudaMalloc((void**)&dev_c, sizeof(float));
	add<<<1,1>>>(1,2,dev_c);  //launch add on GPU
	
	float c;
	cudaMemcpy(&c,dev_c,sizeof(float),cudaMemcpyDeviceToHost);

	cout<<"c = "<<c<<endl;
	return 0;
}
