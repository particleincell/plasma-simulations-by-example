#include <thread>
#include <iostream>
#include <math.h>

using namespace std;
void f(double v, double *res) {(*res) = sqrt(v);}

class B {
public:
  B (double v) : v{v} {};	//save value
  void operator() () {res = sqrt(v);} 
  double res;
protected:
  double v;
};

int main() {
  //thread with a function
  double f_res;  //return value for t1
  thread t1(f,0.2,&f_res); 

  //thread with a reference to a functor class
  B b(2);	
  thread t2(std::ref(b));	

  //thread with a lambda function
  double lambda_res; 
  auto lambda = [&] (double v) {lambda_res = sqrt(v);};
  thread t3(lambda, 4); 
  
  t1.join();	//wait for the threads to finish
  t2.join();
  t3.join();

  cout<<"f_res = "<<f_res<<endl; //result from f
  cout<<"b.res = "<<b.res<<endl; //result from B()
  cout<<"lambda_res = "<<lambda_res<<endl; //result from B()

  return 0;
}
