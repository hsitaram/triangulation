#include<iostream>
#include<fstream>
#include<cmath>
#define N 128
#define PI 3.14159265
int main()
{
	int np;
	double a;
	double b;
	double dth;
	double x,y;
	std::ofstream outfile("polypoints.inp");
	
	a=2.0; b=1.0; 
	np=N;

	dth=2*PI/double(np-1);

	outfile<<np-1<<"\n";

	for(int i=0;i<np-1;i++)
	{
		x=a*cos(i*dth);
		y=b*sin(i*dth);
		outfile<<x<<"\t"<<y<<"\n";	
	}


	outfile.close();

	return(0);
}
