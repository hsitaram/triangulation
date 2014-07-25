#include"polygon.h"
#include"triangulation.h"
#include<cstdlib>

int main(int argc,char *argv[])
{
	std::string option;
	int n;
	double *p;
        polygon pgon;
	triangulation triobj;

	if(argc == 1)
	{
		std::cout<<"\n\n**********************************************\n";
		std::cout<<"usage: ./dlny.exe <algorithm>\n";
		std::cout<<"For centroid insertion use ./dlny.exe cen\n";
		std::cout<<"For Bowyer Watson use ./dlny.exe bw\n";
		std::cout<<"**********************************************\n\n";
	}
	else
	{
		option=argv[1];
		std::ifstream infile;
		infile.open("polypoints.inp");

		infile>>n;

		p = new double[DIM*n];
		for(int i=0;i<n;i++)
		{
			infile>>p[DIM*i]>>p[DIM*i+1];
		}
		infile.close();


		pgon.assignpolypoints(p,n);
		pgon.printpoints();
		pgon.cutear();
		pgon.printtriangles();

		triobj.setpolydomain(&pgon);
		std::cout<<"minlength:"<<pgon.minlength<<"\n";

		if(option=="bw")
		{
			triobj.bwalgorithm(pgon.minlength);
		}
		else if(option=="cen")
		{
			triobj.centroidinsert(pgon.minlength);
		}
		else
		{
			triobj.centroidinsert(pgon.minlength);
		}

		triobj.printtridata();
		triobj.printtrianglesgnuplot();
		triobj.printtrianglesvtu();
	}

	return(0);
}
