#include"edge.h"

bool edge::edgeintersect(edge e)
{
	double x1,y1;
	double x2,y2;
	double x3,y3;
	double x4,y4;

	double px,py;
	double dtr,ntr1,ntr2;
	double dist1,dist2;
	double len1,len2,dotp1,dotp2;

	bool parallel, coincident;

	x1=start[0]; y1=start[1];
	x2=end[0];   y2=end[1];

	e.getendpoints(x3,y3,x4,y4);

	len1 = length;
	len2 = e.getlength();

	dtr = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);

	ntr1 = (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4);
	ntr2 = (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4);

	coincident=false;
	parallel = false;

	if(fabs(dtr) < TOL)
	{
		parallel=true;
		if(fabs(ntr1) < TOL && fabs(ntr2) < TOL)
		{
			coincident=true;
			return(true);
		}	
		else
		{
			return(false);
		}
	}

	px = ntr1/dtr; py=ntr2/dtr;

	dotp1=(px-x1)*(px-x2)+(py-y1)*(py-y2);
	dotp2=(px-x3)*(px-x4)+(py-y3)*(py-y4);

	if(dotp1<=TOL && dotp2<=TOL)
	{
		return(true);
	}
	else
	{
		return(false);
	}
	
}

//==========================================================
bool edge::isitsame(edge e)
{
	bool same;
	double x1,y1,x2,y2;
	double x3,y3,x4,y4;

	double midx1,midy1,midx2,midy2;
	double dtr;
	
	x1=start[0]; y1=start[1];
	x2=end[0];   y2=end[1];
	e.getendpoints(x3,y3,x4,y4);

	dtr = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);

	midx1=0.5*(x1+x2); midx2=0.5*(x3+x4);
	midy1=0.5*(y1+y2); midy2=0.5*(y3+y4);

	same=false;

	if(fabs(dtr) < TOL)
	{
		if(fabs(midx1-midx2) < TOL && fabs(midy1-midy2) < TOL)
		{
			same=true;
		}	
	}

	return(same);
}
//==========================================================
void edge::substitutelrcell(int oldcellnum,int newcellnum)
{
	if(adjcellids[0] == oldcellnum)
	{
		adjcellids[0]=newcellnum;
	}
	else
	{
		adjcellids[1]=newcellnum;
	}
}
//==========================================================
