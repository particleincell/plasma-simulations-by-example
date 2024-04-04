#include <thread>
#include <iostream>
using namespace std;

void add(double a, double b, double *c) {(*c)=a+b;}

int main() 
{
 double c;
 thread t1(add,1,2,&c);
 t1.join();	//wait for the thread to finish
 cout<<"c = "<<c<<endl;
 return 0;
}
