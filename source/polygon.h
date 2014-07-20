#include"triangle.h"
#include"common.h"
class polygon
{
	private:
	
	int numpoints;
	std::vector<int> polyearcutpoints;

	std::vector<int> earpoints;
	std::vector<int> reflexpoints;

	std::vector<edge> polyedges;

	void findearandreflexpoints();

	public:
	
	std::vector<double> allpoints;
	std::vector<triangle> polytriangles;
	double minlength;

	void assignpolypoints(double *p,int n);
	void printpoints();
	void cutear();
	void printtriangles();

};
