#include <thread>
#include <iostream>

using namespace std;
void f(double v) {
  cout<<"Hello from function f("<<v<<")"<<endl;
}

class B {
public:
  B(double v): v{v} {}
  void operator() () {
    cout<<"Greeting from class B with v="<<v<<endl;} 
  double v;
};

int main() {
  thread t1(f,0.1); //create 3 threads
  thread t2(B(0.2));
  B b(0.0);
  b.v = 0.3;	//reset value
  thread t3(b);	//another approach
  
  t1.join();	//wait for the threads to finish
  t2.join();
  t3.join();

  return 0;
}
