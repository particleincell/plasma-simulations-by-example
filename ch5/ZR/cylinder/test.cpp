#include"World.h"
#include<iostream>
#include<math.h>
#include<vector>
using namespace std;
Rnd rnd;

int main()
{
  vector<int> hist(10);
for (int i=0;i<100000;i++) 
{
  double r=sqrt(rnd());
  int b = (int)(r*10);
  hist[b]++;
}
	
for (int b=0;b<10;b++)
 {
	for (int i=0;i<hist[b]/1000;i++)
	 cout<<"*";
      cout<<endl;	
 }
	return 0;
}
