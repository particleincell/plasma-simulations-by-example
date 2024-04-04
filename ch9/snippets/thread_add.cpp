#include <thread>
#include <iostream>
using namespace std;

int main() {
	for (int k=1;k<=10000;k++)  //run 10,000 times
	{
		int i=0;
		auto f = [&] {++i;};  //lambda expression to increment i
		thread t1(f); 		  //create two threads
		thread t2(f);
		t1.join(); t2.join(); //wait for threads to finish
		cout<<(i==2?'.':'*'); //output . if i=2, * otherwise
		if (k%100==0) cout<<endl;
	}
	return 0;
}

