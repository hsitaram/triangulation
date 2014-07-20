#include"triangle.h"

class triangulation
{
	int numnodes;
	std::vector<double> nodes;
	std::vector<edge> edgelist;
	std::vector<triangle> trianglelist;
	
	bool performflips();
	bool isedgepresent(edge e,int &pos);
	void addedgetolist(edge e,int trinum);
	void reordertriedges();
	void addnewpoint(int tid);
	public:

	void assign_nodes(std::vector<double> points){ nodes=points; numnodes=nodes.size()/2;}
	void assign_triangles(std::vector<triangle> alltriangles);

	void printtridata();
	void printtrianglesgnuplot();
	void printtrianglesvtu();
	void delaunay(double minsidelen);

};
