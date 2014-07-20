#include"polygon.h"
#include"triangulation.h"

int main()
{
	int n;
	double *p;
        polygon pgon;
	triangulation triobj;

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

	triobj.assign_nodes(pgon.allpoints);
	triobj.assign_triangles(pgon.polytriangles);

	std::cout<<"minlength:"<<pgon.minlength<<"\n";
	triobj.delaunay(pgon.minlength);

	triobj.printtridata();
	triobj.printtrianglesgnuplot();
	triobj.printtrianglesvtu();


	return(0);
}
