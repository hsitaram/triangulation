#include"triangulation.h"
//=======================================================================================
void triangulation::assign_triangles(std::vector<triangle> alltriangles)
{
	trianglelist=alltriangles;
	formedgelist();


}
//=======================================================================================
void triangulation::formedgelist()
{
	int tp[3];
	edge t1,t2,t3;

	edgelist.clear();	
	edgelist.resize(0);

	//resize edge ids because we use push_back later
	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		trianglelist[i].edgepos=0;
	}

	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		trianglelist[i].getedges(t1,t2,t3);
		trianglelist[i].getpointids(tp);

		addedgetolist(t1,i);
		addedgetolist(t2,i);	
		addedgetolist(t3,i);
	}


}
//=======================================================================================
bool triangulation::isedgepresent(edge e,int &pos)
{
	bool present;

	present=false;
	pos=-1;

	for(unsigned int i=0;i<edgelist.size();i++)
	{
		if(edgelist[i].isitsame(e))
		{
			present=true;
			pos=i;
			break;
		}
	}

	return(present);
}
//=======================================================================================
void triangulation::addedgetolist(edge e,int trinum)
{
	int pos;
	bool present;
	int edgenum;

	present=isedgepresent(e,pos);

	if(!present)
	{
		e.adjcellids[0]=trinum;
		edgelist.push_back(e);
		edgenum=edgelist.size()-1;
		trianglelist[trinum].edgeids[trianglelist[trinum].edgepos++]=edgenum;

	}
	else
	{
		edgelist[pos].adjcellids[1]=trinum;
		trianglelist[trinum].edgeids[trianglelist[trinum].edgepos++]=pos;
	}

}
//=======================================================================================
void triangulation::centroidinsert(double minsidelen)
{
	int it,itmax,i;
	bool flipflag;
	int maxareatri;
	double maxarea;
	double avgsidelen;
	
	itmax=50;

	for(it=0;it<itmax;it++)
	{
		flipflag=performflips();
		if(!flipflag)
		{
			break;
		}

	}

	for(i=0;i<itmax*100;i++)
	{
		reordertriedges();

		maxarea=0.0;
		avgsidelen=0.0;
		maxareatri=0;
		for(unsigned int j=0;j<trianglelist.size();j++)
		{
			avgsidelen=avgsidelen+trianglelist[j].avglen;
			if(trianglelist[j].area > maxarea)
			{
				maxarea=trianglelist[j].area;
				maxareatri=j;
			}	
		}

		avgsidelen=avgsidelen/double(trianglelist.size());

		if(avgsidelen <= minsidelen)
		{
			break;
		}

		addnewpoint(maxareatri);

		for(it=0;it<itmax;it++)
		{
			flipflag=performflips();
			if(!flipflag)
			{
				break;
			}

		}
		//std::cout<<"flip count:"<<it<<"\n";
	}

	std::cout<<"avgsidelen:"<<avgsidelen<<"\n";
	std::cout<<"number of triangles:"<<trianglelist.size()<<"\n";
	std::cout<<"number of nodes:"<<nodes.size()/2<<"\n";
	std::cout<<"number of edges:"<<edgelist.size()<<"\n";
	std::cout<<"no: of iterations:"<<i<<"\n";

}
//=======================================================================================
void triangulation::printtridata()
{
	int *p;
	std::ofstream outfile("tri.msh");

	//nodes
	for(int i=0;i<numnodes;i++)
	{
		outfile<<i<<"\t"<<nodes[2*i]<<"\t"<<nodes[2*i+1]<<"\n";
	}

	outfile<<"\n\n";

	//triangles
	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		p = trianglelist[i].getpointids();
		outfile<<i<<"\t"<<p[0]<<"\t"<<p[1]<<"\t"<<p[2]<<"\n";
	}
	outfile<<"\n\n";

	//edges
	for(unsigned int i=0;i<edgelist.size();i++)
	{
		outfile<<i<<"\t"<<edgelist[i].nodeids[0]<<"\t"<<edgelist[i].nodeids[1]<<"\t"
			<<edgelist[i].adjcellids[0]<<"\t"<<edgelist[i].adjcellids[1]<<"\n";
	}

	outfile<<"\n\n";	

	//triangle edges
	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		outfile<<i<<"\t"<<trianglelist[i].edgeids[0]<<"\t"
			<<trianglelist[i].edgeids[1]
			<<"\t"<<trianglelist[i].edgeids[2]<<"\n";
	}
	outfile<<"\n\n";

	outfile.close();

}
//=======================================================================================
void triangulation::printtrianglesgnuplot()
{
	int ind1,ind2;
	std::ofstream outfile("tri.dat");

	for(unsigned int i=0;i<edgelist.size();i++)
	{
		ind1=edgelist[i].nodeids[0];
		ind2=edgelist[i].nodeids[1];
		outfile<<nodes[2*ind1]<<"\t"<<nodes[2*ind1+1]<<"\n";
		outfile<<nodes[2*ind2]<<"\t"<<nodes[2*ind2+1]<<"\n\n\n";
	}

	outfile.close();

}
//=======================================================================================
bool triangulation::performflips()
{
	triangle t1,t2;
	edge e,newedge;
	int lcell,rcell;
	bool boundaryedgeflag;
	int cl,cr; //corners
	double xcl,ycl,xcr,ycr;
	int nd1,nd2;
	int otheredgesL[2],otheredgesR[2];

	int nedges;
	bool flipflag;

	nedges=edgelist.size();
	//nedges=16;

	flipflag=false;
	for(int i=0;i<nedges;i++)
	{
		e=edgelist[i];


		nd1=e.nodeids[0]; nd2=e.nodeids[1];

		lcell=e.adjcellids[0];
		rcell=e.adjcellids[1];
		

		boundaryedgeflag=((lcell==-1) || (rcell==-1))?true:false;

		if(!boundaryedgeflag)
		{
			cl=trianglelist[lcell].getoppositepoint(e);
			cr=trianglelist[rcell].getoppositepoint(e);

			xcl = nodes[2*cl]; ycl=nodes[2*cl+1];
			xcr = nodes[2*cr]; ycr=nodes[2*cr+1];

			if(trianglelist[lcell].ispointinsidecircumcircle(xcr,ycr) || 
					trianglelist[rcell].ispointinsidecircumcircle(xcl,ycl))
			{
				flipflag=true;

				//get other edges
				trianglelist[lcell].getotheredgeids(i,otheredgesL);
				trianglelist[rcell].getotheredgeids(i,otheredgesR);

				//create new triangles
				t1.setvertices(nodes[2*cl],nodes[2*cl+1],
						nodes[2*nd1],nodes[2*nd1+1],
						nodes[2*cr],nodes[2*cr+1],cl,nd1,cr,"setting t1");

				t2.setvertices(nodes[2*cr],nodes[2*cr+1],
						nodes[2*nd2],nodes[2*nd2+1],
						nodes[2*cl],nodes[2*cl+1],cr,nd2,cl,"setting t2");

				trianglelist[lcell] = t1;
				trianglelist[rcell] = t2;

				//create new edge
				newedge.setedgeparams(nodes[2*cl],nodes[2*cl+1],
						nodes[2*cr],nodes[2*cr+1],cl,cr);
				newedge.adjcellids[0]=lcell;
				newedge.adjcellids[1]=rcell;

				edgelist[i] = newedge;

				//add edgeids to new triangles
				trianglelist[lcell].edgeids[trianglelist[lcell].edgepos++]=i;
				trianglelist[rcell].edgeids[trianglelist[rcell].edgepos++]=i;

				//take care of surrounding edges
				if(edgelist[otheredgesL[0]].nodeids[0] == nd2 || edgelist[otheredgesL[0]].nodeids[1] == nd2)
				{
					edgelist[otheredgesL[0]].substitutelrcell(lcell,rcell);
					trianglelist[rcell].edgeids[trianglelist[rcell].edgepos++]=otheredgesL[0];

					edgelist[otheredgesL[1]].substitutelrcell(lcell,lcell);
					trianglelist[lcell].edgeids[trianglelist[lcell].edgepos++]=otheredgesL[1];
				}
				else
				{
					edgelist[otheredgesL[0]].substitutelrcell(lcell,lcell);
					trianglelist[lcell].edgeids[trianglelist[lcell].edgepos++]=otheredgesL[0];

					edgelist[otheredgesL[1]].substitutelrcell(lcell,rcell);
					trianglelist[rcell].edgeids[trianglelist[rcell].edgepos++]=otheredgesL[1];

				}

				if(edgelist[otheredgesR[0]].nodeids[0] == nd2 || edgelist[otheredgesR[0]].nodeids[1] == nd2)
				{
					edgelist[otheredgesR[0]].substitutelrcell(rcell,rcell);
					trianglelist[rcell].edgeids[trianglelist[rcell].edgepos++]=otheredgesR[0];

					edgelist[otheredgesR[1]].substitutelrcell(rcell,lcell);
					trianglelist[lcell].edgeids[trianglelist[lcell].edgepos++]=otheredgesR[1];
				}
				else
				{
					edgelist[otheredgesR[0]].substitutelrcell(rcell,lcell);
					trianglelist[lcell].edgeids[trianglelist[lcell].edgepos++]=otheredgesR[0];

					edgelist[otheredgesR[1]].substitutelrcell(rcell,rcell);
					trianglelist[rcell].edgeids[trianglelist[rcell].edgepos++]=otheredgesR[1];
				}


			}	
		}

	}

	return(flipflag);	

}
//=======================================================================================
void triangulation::printtrianglesvtu()
{
	char underscore;
	int offset;
	double zero=0.0;
	double dval;
	int intval;
	int p[3];

	std::ofstream outfile;
	outfile.open("triangles.vtu");

	offset=0;
	underscore='_';

	outfile<<"<?xml version=\"1.0\"?>\n";	
	outfile<<"<VTKFile type=\"UnstructuredGrid\" ";
	outfile<<"version=\"0.1\" byte_order=\"LittleEndian\">\n";	
	outfile<<"<UnstructuredGrid>\n";
	outfile<<"<Piece NumberOfPoints=\""<<numnodes<<"\" NumberOfCells=\""<<trianglelist.size()<<"\">\n";

	outfile<<"<CellData>\n";
	outfile<<"<DataArray type=\"Float64\" Name=\"circumrad\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
	outfile<<"</CellData>\n";

	offset = offset + sizeof(int) + sizeof(double)*trianglelist.size();
	outfile<<"<Points>\n";
	outfile<<"<DataArray name=\"Position\" type=\"Float64\" ";
	outfile<<"NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
	outfile<<"</Points>\n";

	offset = offset + sizeof(int) + 3*sizeof(double)*numnodes;

	outfile<<"<Cells>\n";
	outfile<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset<<"\"/>\n";

	offset = offset + sizeof(int) + 3*sizeof(int)*trianglelist.size();
	outfile<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset<<"\"/>\n";

	offset = offset + sizeof(int) + sizeof(int)*trianglelist.size();	
	outfile<<"<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\""<<offset<<"\"/>\n";

	outfile<<"</Cells>\n";
	outfile<<"</Piece>\n";
	outfile<<"</UnstructuredGrid>\n";
	outfile<<"<AppendedData encoding=\"raw\">\n";

	outfile.write(&underscore,sizeof(char));

	intval=sizeof(double)*trianglelist.size();
	outfile.write((char *)&intval,sizeof(int));
	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		outfile.write((char *)&trianglelist[i].circrad,sizeof(double));	
	}

	intval=3*sizeof(double)*numnodes;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<numnodes;i++)
	{
		dval=nodes[2*i];
		outfile.write((char *)&dval,sizeof(double));
		dval=nodes[2*i+1];
		outfile.write((char *)&dval,sizeof(double));
		outfile.write((char *)&zero,sizeof(double));
	}

	intval=3*sizeof(int)*trianglelist.size();
	outfile.write((char *)&intval,sizeof(int));
	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		trianglelist[i].getpointids(p);
		outfile.write((char *)&p[0],sizeof(int));	
		outfile.write((char *)&p[1],sizeof(int));	
		outfile.write((char *)&p[2],sizeof(int));	
	}

	intval=sizeof(int)*trianglelist.size();
	outfile.write((char *)&intval,sizeof(int));
	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		intval=(i+1)*3;
		outfile.write((char *)&intval,sizeof(int));
	}	
	intval=sizeof(int)*trianglelist.size();
	outfile.write((char *)&intval,sizeof(int));
	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		intval=5;
		outfile.write((char *)&intval,sizeof(int));
	}	

	outfile.close();

}
//=======================================================================================
void triangulation::reordertriedges()
{
	int id1,id2;
	int pids[3],k;
	triangle t;
	edge test_e;
	std::vector<int> edgeids;

	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		t=trianglelist[i];
		t.getpointids(pids);

		edgeids.resize(0);

		for(int d=0;d<3;d++)
		{

			id1=pids[d];
			id2=pids[(d+1)%3];

			test_e.setedgeparams(nodes[2*id1],nodes[2*id1+1],
					nodes[2*id2],nodes[2*id2+1],id1,id2);

			for(k=0;k<3;k++)
			{
				if(edgelist[t.edgeids[k]].isitsame(test_e))
					break;
			}

			edgeids.push_back(t.edgeids[k]);

		}

		for(int k=0;k<3;k++)
		{
			trianglelist[i].edgeids[k]=edgeids[k];
		}

	}

	edgeids.clear();
}
//=======================================================================================
void triangulation::addnewpoint(int tid)
{
	edge e[3]; //new edges
	int te[3]; // tri edges
	int tp[3]; // tri points
	int tp1,tp2; //tri points
	int nnid; //new node id
	int newtid[3];
	int neweid[3];

	triangle t[3];
	triangle told;

	trianglelist[tid].getpointids(tp);
	told = trianglelist[tid];

	te[0]=trianglelist[tid].edgeids[0];
	te[1]=trianglelist[tid].edgeids[1];
	te[2]=trianglelist[tid].edgeids[2];

	//add point
	nodes.push_back(trianglelist[tid].centroid[0]);
	nodes.push_back(trianglelist[tid].centroid[1]);
	numnodes=numnodes+1;
	nnid = numnodes-1;

	//create new triangles
	for(int k=0;k<3;k++)
	{
		tp1=tp[k];
		tp2=tp[(k+1)%3];

		t[k].setvertices(nodes[2*tp1],nodes[2*tp1+1],
				nodes[2*nnid],nodes[2*nnid+1],
				nodes[2*tp2],nodes[2*tp2+1],tp1,nnid,tp2);

	}

	trianglelist[tid]=t[0];
	newtid[0]=tid;

	trianglelist.push_back(t[1]);
	trianglelist.push_back(t[2]);

	newtid[1]=trianglelist.size()-2;
	newtid[2]=trianglelist.size()-1;

	//add edges
	for(int k=0;k<3;k++)
	{
		tp1=tp[k];

		//std::cout<<"tp1:"<<tp1<<"\n";

		e[k].setedgeparams(nodes[2*tp1],nodes[2*tp1+1],
				nodes[2*nnid],nodes[2*nnid+1],tp1,nnid);

		e[k].adjcellids[0]=newtid[k]; 
		e[k].adjcellids[1]=newtid[(k+5)%3];

	}

	edgelist.push_back(e[0]);
	edgelist.push_back(e[1]);
	edgelist.push_back(e[2]);

	neweid[0]=edgelist.size()-3; 
	neweid[1]=edgelist.size()-2;
	neweid[2]=edgelist.size()-1;

	//correct r and l cell	
	edgelist[te[0]].substitutelrcell(tid,newtid[0]);
	edgelist[te[1]].substitutelrcell(tid,newtid[1]);
	edgelist[te[2]].substitutelrcell(tid,newtid[2]);

	//update edges for new triangles

	for(int k=0;k<3;k++)
	{
		trianglelist[newtid[k]].edgeids[trianglelist[newtid[k]].edgepos++] = neweid[k];
		trianglelist[newtid[k]].edgeids[trianglelist[newtid[k]].edgepos++] = neweid[(k+1)%3];
		trianglelist[newtid[k]].edgeids[trianglelist[newtid[k]].edgepos++] = te[k];
	}	



}
//=======================================================================================
bool triangulation::insertpointbowyerwatson(int tid)
{
	int ntid; //neighboring triangle id
	int neighbortri[3];
	bool neighbordeletionflaglist[3];
	int tp[3];
	int te[3];
	double newpoint[2];
	int ntrinodes[3];
	int nnid,tp1,tp2;
	std::vector<int> polynodes;
	std::vector<int> deltrilist;
	std::vector<triangle> addtrilist;
	int npolynodes;
	bool insertflag=false;

	triangle t;

	deltrilist.resize(0);
	addtrilist.resize(0);
	polynodes.resize(6);
	trianglelist[tid].getpointids(tp);

	neighbordeletionflaglist[0]=false;
	neighbordeletionflaglist[1]=false;
	neighbordeletionflaglist[2]=false;

	te[0]=trianglelist[tid].edgeids[0];
	te[1]=trianglelist[tid].edgeids[1];
	te[2]=trianglelist[tid].edgeids[2];

	ntrinodes[0]=-1; 
	ntrinodes[1]=-1; 
	ntrinodes[2]=-1;
	

	newpoint[0]=trianglelist[tid].circ_x; 
	newpoint[1]=trianglelist[tid].circ_y;

	if(polydomain->ispointinside(newpoint[0],newpoint[1]))
	{
		for(int i=0;i<3;i++)
		{
			neighbortri[i]=(edgelist[te[i]].adjcellids[0]==tid)?
				edgelist[te[i]].adjcellids[1]:edgelist[te[i]].adjcellids[0];

			ntid=neighbortri[i];

			if(ntid != -1)
			{	
				ntrinodes[i]=trianglelist[ntid].getoppositepoint(edgelist[te[i]]);
			}
		}

		for(int i=0;i<3;i++)
		{
			ntid=neighbortri[i];
			if(ntid != -1)
			{
				if(trianglelist[ntid].ispointinsidecircumcircle(newpoint[0],newpoint[1]))
				{
					neighbordeletionflaglist[i]=true;
				}
			}
		}

		polynodes[0] = tp[0];
		polynodes[1] = ntrinodes[0];
		polynodes[2] = tp[1];
		polynodes[3] = ntrinodes[1];
		polynodes[4] = tp[2];
		polynodes[5] = ntrinodes[2];

		//delete from behind to leave indexing unaffected.	
		for(int i=2;i>=0;i--)
		{
			if(!neighbordeletionflaglist[i])
			{
				polynodes.erase(polynodes.begin()+2*i+1);
			}
		}
		
		nnid=numnodes;
		//std::cout<<"ntri points:"<<ntrinodes[0]<<"\t"<<ntrinodes[1]<<"\t"<<ntrinodes[2]<<"\n";

		//delete triangles
		deltrilist.push_back(tid);
		for(int i=0;i<3;i++)
		{
			if(neighbordeletionflaglist[i])
			{
				ntid=neighbortri[i];
				deltrilist.push_back(ntid);
			}
		}

		//form triangles and add to list
		npolynodes=polynodes.size();

		//std::cout<<"npolynodes:"<<npolynodes<<"\n";
		addtrilist.resize(0);
		insertflag=true;
		for(unsigned int i=0;i<polynodes.size();i++)
		{
			tp1 = polynodes[i];
			tp2 = polynodes[(i+1)%npolynodes];

			t.setvertices(nodes[2*tp1],nodes[2*tp1+1],
					newpoint[0],newpoint[1],
					nodes[2*tp2],nodes[2*tp2+1],tp1,nnid,tp2);

			if(t.area < TOL)
			{
				insertflag=false;
				break;
			}
			else
			{
				addtrilist.push_back(t);
			}
		}

		if(insertflag)
		{
			deletetrifromlist(deltrilist);
			for(unsigned int i=0;i<addtrilist.size();i++)
			{
				trianglelist.push_back(addtrilist[i]);
			}

			nodes.push_back(newpoint[0]);
			nodes.push_back(newpoint[1]);
			numnodes=numnodes+1;
		}

	}

	//local cleanup
	addtrilist.clear();
	polynodes.clear();
	deltrilist.clear();
	return(insertflag);

}
//===================================================================================================
void triangulation::bwalgorithm(double minsidelen)
{
	int it,itmax=10;
	bool flipflag;
	int maxareatri;
	double maxarea;
	double avgsidelen;
	bool insertflag;
	bool alltriangleschecked;
	std::vector<int> marked;

	for(it=0;it<itmax;it++)
	{
		flipflag=performflips();
		if(!flipflag)
		{
			break;
		}

	}

	for(int i=0;i<itmax*300;i++)
	{
		reordertriedges();
		insertflag=false;
		marked.resize(0);
		alltriangleschecked=false;

		while(!insertflag)
		{
			maxarea=0.0;
			maxareatri=0;

			if(marked.size() > trianglelist.size()-1)
			{
				alltriangleschecked=true;
				break;
			}

			for(unsigned int j=0;j<trianglelist.size();j++)
			{
				if(trianglelist[j].area > maxarea && 
						(std::find(marked.begin(),marked.end(),j) == marked.end()))
				{
					maxarea=trianglelist[j].area;
					maxareatri=j;
				}	
			}

			marked.push_back(maxareatri);
			insertflag = insertpointbowyerwatson(maxareatri);
		}

		if(alltriangleschecked)
		{
			break;
		}
		//testformedges();
		formedgelist();
		//reordertriedges();

		for(it=0;it<itmax;it++)
		{
			flipflag=performflips();
			if(!flipflag)
			{
				break;
			}

		}
		avgsidelen=0.0;
		for(unsigned int j=0;j<trianglelist.size();j++)
		{
			avgsidelen=avgsidelen+trianglelist[j].avglen;
		}
		avgsidelen = avgsidelen/double(trianglelist.size());

		if(avgsidelen <= minsidelen)
		{
			break;
		}

		marked.clear();
		std::cout<<"flip count:"<<it<<"\titer:"<<i<<"\tavgsidelen:"<<avgsidelen<<
			"\tminsidelen:"<<minsidelen<<"\n";
	}
	
	std::cout<<"number of triangles:"<<trianglelist.size()<<"\n";
	std::cout<<"number of nodes:"<<nodes.size()/2<<"\n";
	std::cout<<"number of edges:"<<edgelist.size()<<"\n";
	
	//cleanup
	marked.clear();

}
//===================================================================================================
void triangulation::deletetrifromlist(std::vector<int> deltrilist)
{
	int numtri;
	std::vector<triangle> newtrilist;

	numtri=trianglelist.size();

	newtrilist.resize(0);

	for(int i=0;i<numtri;i++)
	{
		if(std::find(deltrilist.begin(),deltrilist.end(),i)==deltrilist.end())
		{
			newtrilist.push_back(trianglelist[i]);
		}
	}

	trianglelist.clear();
	trianglelist=newtrilist;

	//local clean up
	newtrilist.clear();

}
//===================================================================================================
void triangulation::testformedges()
{
	for(int i=0;i<numnodes;i++)
	{
		std::cout<<i<<"\t"<<nodes[2*i]<<"\t"<<nodes[2*i+1]<<"\n";
	}


	std::cout<<"printing tri.msh\n";
	printtridata();

	formedgelist();
	reordertriedges();
	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		std::cout<<i<<"\t"<<trianglelist[i].getpointids()[0]<<"\t"<<
		trianglelist[i].getpointids()[1]<<"\t"<<trianglelist[i].getpointids()[2]<<"\n";
	}

}
//===================================================================================================
