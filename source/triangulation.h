#include"triangle.h"
#include"polygon.h"
#include<algorithm>

class triangulation
{
	int numnodes;
	std::vector<double> nodes;
	std::vector<edge> edgelist;
	std::vector<triangle> trianglelist;
	polygon *polydomain;
	
	bool performflips();
	bool isedgepresent(edge e,int &pos);
	void addedgetolist(edge e,int trinum);
	void reordertriedges();
	void addnewpoint(int tid);
	bool insertpointbowyerwatson(int tid);
	void formedgelist();
	void deletetrifromlist(std::vector<int> deltrilist);

	void testformedges();	
	void assign_nodes(std::vector<double> points)
	{ 
		nodes=points; 
		numnodes=nodes.size()/2;
	}
	void assign_triangles(std::vector<triangle> alltriangles);
	
	public:

        void setpolydomain(polygon *p)
	{
		polydomain=p;
		assign_nodes(polydomain->allpoints);
		assign_triangles(polydomain->polytriangles);
	}	

	void printtridata();
	void printtrianglesgnuplot();
	void printtrianglesvtu();
	void centroidinsert(double minsidelen);
	void bwalgorithm(double minsidelen);

};
