
#pragma once

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>

using std::vector;
typedef double EdgeCost;
//THE LINE BELOW WAS SET BY THOMAS
#define EDGE_COST_PRECISION 0.0000001 
//#define EDGE_COST_PRECISION 0.0001
#define INFINITE_COST 1e32
#define DEBUG_VERTEX -1

int oldindex[10000005];
int oldedgeid[10000005];
int newindex[10000005];


struct EdgeDescriptor {
	int v;
	int w;
	EdgeCost cost;

	EdgeDescriptor (int _v, int _w, EdgeCost _cost) {
		if (_v < _w) {
			v = _v;
			w = _w;
		} else {
			v = _w;
			w = _v;
		}
		cost = _cost;
	}
};


struct GraphDescriptor {
	int m;
	int n;
	int tcount;

	vector<EdgeDescriptor> edges;
	vector<bool> terminal;
	vector<int> Belong_to_Patition;
	double fixed;

	GraphDescriptor () {
		// add dummy edge at position 0; edge IDs start at 1
		edges.clear();
		edges.push_back(EdgeDescriptor(-1,-1,-1));
		n = m = -1;
		tcount = 0;
		fixed=0;
	}
	inline void SetFixed(double f) {fixed = f;}
	inline void IncFixed(double f) {fixed += f;}

	inline int GetEdgeCount() const {
		return (int)(edges.size() - 1);//m;
	}

	inline void CheckRange (int v) {
		if (v<=0 || v>n) {
			fprintf (stderr, "Vertex %d is out of range [1,%d].\n", v, n);
			exit(-1);
		}
	}

	void SetVertices(int _n) {
        n = _n; 
        terminal.resize(n+1,false); 
        Belong_to_Patition.resize(n+1,0);
		/*for(int i=1;i<=n;i++)
		{
			terminal[i]=false;
			Belong_to_Patition[i]=0;
		}*/
    }

	void SetEdges(int _m) {
	}

	int AddEdge(int v, int w, EdgeCost cost) {
		int cursize = (int)edges.size();
		CheckRange(v);
		CheckRange(w);
		edges.push_back(EdgeDescriptor(v,w,cost));
		return cursize;
	}

	void MakeTerminal(int v) {
        CheckRange(v); 
        if (!terminal[v]) {
            terminal[v]=true; 
            tcount++;
        }
    }

	void UnmakeTerminal(int v) {
		CheckRange(v); 
        if (terminal[v]) {
            terminal[v] = false; 
            tcount--;
        }
	}

	int TerminalCount() {
        return tcount;
    }

	void SetBelongPartition(int v,int flag){
		Belong_to_Patition[v]=flag;
	}
	int GetBelong_to_Patition(int v){
		return Belong_to_Patition[v];
	}
};


struct SPGArc {
	int label; //arc label
	int head;
	EdgeCost cost;

	inline void Set (int h, EdgeCost c, int l) {
        head = h; 
        cost = c; 
        label = l;
    }
};


// dynamic graph, with new edges added as they are created
struct DynamicGraph {
	int m, n, tcount;
	vector <SPGArc> arclist;
	vector <int> first;
	vector <int> last;
	GraphDescriptor gd;

	void OutputAdjacencyList (FILE *file, int v) {
		fprintf (stderr, "%4d:", v);
		SPGArc *a, *end;
		for (GetBounds(v,a,end); a<end; a++) {
			fprintf (stderr, " %d:%.0f:%d ", a->head, a->cost, a->label);
		}
		fprintf (stderr, "\n");
	}
	
	void OutputGraph (FILE *file) {
		for (int v=1; v<=n; v++) 
            OutputAdjacencyList(file, v);
	} 

	inline double GetFixedCost() {return gd.fixed;}
	inline void SetFixedCost(EdgeCost f) {gd.SetFixed(f);}
	inline void IncFixedCost(EdgeCost f) {gd.IncFixed(f);}

	inline EdgeCost GetCost (int e) const {return gd.edges[e].cost;}
	inline int VertexCount() const {return n;}
	inline int EdgeCount() const {return gd.GetEdgeCount();}
	inline int TerminalCount() const {return gd.tcount;}
	inline int GetDegree(int v) const {return last[v]-first[v];}
	inline int GetReservedDegree(int v) const {return first[v+1] - first[v];}
	inline int GetFirstEdgeLabel(int v) const { 
		return arclist[first[v]].label;
	}

	inline bool IsTerminal(int v) const {return gd.terminal[v];}

	inline void GetBounds (int v, SPGArc * &start, SPGArc * &end) {
		start = &arclist[first[v]];
		end = &arclist[last[v]];
	}
	inline void SetVertices(int _n) {gd.SetVertices(_n);}
	inline void SetEdges(int _m) {gd.SetEdges(_m);} 
	inline int BatchAddEdge(int v, int w, EdgeCost cost) {return gd.AddEdge(v,w,cost);}
	inline void MakeTerminal(int v) {gd.MakeTerminal(v);}
	inline void UnmakeTerminal(int v) {gd.UnmakeTerminal(v);}
	inline void GetEndpoints (int e, int &v, int &w) {
		v = gd.edges[e].v;
		w = gd.edges[e].w;
	}

	inline int InsertEdge (int _v, int _w, EdgeCost cost) {
		int newid = gd.AddEdge(_v,_w,cost);
		int endpoint[2] = {_v,_w};

		// insert new arc at the back of each adjacency list
		for (int i=0; i<2; i++) {
			int v = endpoint[i];
			int w = endpoint[1-i];
			arclist[last[v]].Set(w,cost,newid);
			last[v]++;
		}
		return newid;
	}

	inline void HideEdge (int e) {
		const bool verbose = false;
		int neighbor[2];
		GetEndpoints(e,neighbor[0],neighbor[1]);
		
		for (int i=0; i<2; i++) {
			int v = neighbor[i];
			int p;
			for (p = first[v]; p<last[v]; p++) {
				if (arclist[p].label == e) break;
			}

			last[v] --;
			for ( ; p<last[v]; p++) {arclist[p] = arclist[p+1];}
		}
	}

	//--------------------------------------------------
	// Create adjacency lists from plain lists of edges
	//--------------------------------------------------
	void Commit () {
		n = gd.n;
		int ecount = (int)gd.edges.size() - 1; //we start at 1...
		m = gd.GetEdgeCount();

		tcount = gd.tcount;
		arclist.resize(2*m+1);
		first.resize(n+2,0);
		last.resize(n+1,0);

		int v, e;
		// now first will correspond to the degree
		for (e=1; e<=m; e++) {
			first[gd.edges[e].v]++;
			first[gd.edges[e].w]++;
		}

		// convert degrees into first positions
		int acc = 0; //accumulated degree up to this point
		for (v=1; v<=n; v++) {
			int temp = first[v];
			first[v] = acc;
			acc += temp;
		}

		// insert the actual arcs
		for (e=1; e<=m; e++) {
			int v = gd.edges[e].v;
			int w = gd.edges[e].w;
			EdgeCost cost = gd.edges[e].cost;
			arclist[first[v]++].Set(w,cost,e);
			arclist[first[w]++].Set(v,cost,e);
		}

		// restore first pointers (including sentinel)
		for (v = n+1; v > 0; v--) {first[v] = first[v-1];}

		// saved last for now
		for (v=1; v<=n; v++) {last[v] = first[v+1];}

	}

};

struct Graph {
	int m, n, tcount;
	vector <SPGArc> arclist;
	vector <int> first;
	GraphDescriptor gd;
	//vector<int>oldindex;
	//vector<int>newindex;

	Graph () {
		n = m = tcount = 0;
	}

	inline double GetFixedCost() {return gd.fixed;}
	inline void SetFixedCost(double f) {gd.SetFixed(f);}
	inline void IncFixedCost(double f) {gd.IncFixed(f);}

	inline void SetVertices(int _n) {gd.SetVertices(_n);}
	inline void SetEdges(int _m) {gd.SetEdges(_m);}
	inline int AddEdge(int v, int w, EdgeCost cost) {return gd.AddEdge(v,w,cost);}
	inline void MakeTerminal(int v) {gd.MakeTerminal(v);}
	inline void UnmakeTerminal(int v) {gd.UnmakeTerminal(v);}

	inline void SetBelongPartition(int v,int flag) {gd.SetBelongPartition(v,flag);}
	inline int GetBelong_to_Patition(int v) {return gd.GetBelong_to_Patition(v);};

	//--------------------------------------------------
	// Create adjacency lists from plain lists of edges
	//--------------------------------------------------
	void Commit () {
		n = gd.n;
		int ecount = (int)gd.edges.size() - 1;
		m = gd.GetEdgeCount();

		tcount = gd.tcount;
		arclist.resize(2*m+1); 
		first.resize(n+2,0);

		int v, e;
		for (e=1; e<=m; e++) {
			first[gd.edges[e].v]++;
			first[gd.edges[e].w]++;
		}

		int acc = 0;
		for (v=1; v<=n; v++) {
			int temp = first[v];
			first[v] = acc;
			acc += temp;
		}
		for (e=1; e<=m; e++) {
			int v = gd.edges[e].v;
			int w = gd.edges[e].w;
			EdgeCost cost = gd.edges[e].cost;
			arclist[first[v]++].Set(w,cost,e);
			arclist[first[w]++].Set(v,cost,e);
		}
		for (v = n+1; v > 0; v--) {first[v] = first[v-1];}
	}

	inline int GetArcTail (int a) const {
		return (a<=m) ? GetFirstEndpoint(a) : GetSecondEndpoint(a-m);
	}

	inline void GetArcEndpoints(int a, int &v, int &w) {
		if (a <= m) {
			GetEndpoints(a,v,w); //[1..m]: same as undirected edge (smaller first)
		} else {
			GetEndpoints(a-m,w,v); //high: oppositive as undirected edge (higher first)
		}
	}

    inline int GetOutgoingLabel (int v, SPGArc *a) const {
		int w = a->head;
		return (v<w) ? a->label : a->label + m;
	}

	inline int GetIncomingLabel (int v, SPGArc *a) const {
		int w = a->head;
		return (v<w) ? a->label+m : a->label;
	}

    inline int GetIncomingLabel (int v, int w, int alabel) const {
		return (v<w) ? alabel + m : alabel;
	}



	inline int VertexCount() const {return n;}
	inline int EdgeCount() const {return m;}
	inline int TerminalCount() const {return gd.tcount;}
	inline int GetDegree(int v) const {return first[v+1]-first[v];}

	inline bool IsTerminal(int v) const {return gd.terminal[v];}

	inline void GetBounds (int v, SPGArc * &start, SPGArc * &end) {
		start = &arclist[first[v]];
		end = &arclist[first[v+1]];
	}

	// copy current edge costs to "costs"
	void RetrieveCosts (vector<EdgeCost> &costs) {
		costs.clear();
		costs.resize(m+1);
		int esize = (int)gd.edges.size();
		for (int e=0; e<=m; e++) {
			costs[e] = gd.edges[e].cost;
		}
	}

	// update all costs in the graph according 
	void ApplyCosts (vector<EdgeCost> &newcost) {
		for (int v=1; v<=n; v++) {
			SPGArc *a, *end;
			for (GetBounds(v,a,end); a<end; a++) {
				a->cost = newcost[a->label];
			}
		}
		for (int e=0; e<=m; e++) {
			gd.edges[e].cost = newcost[e];
		}
	}

	void OutputGraph (FILE *file) {
		for (int v=1; v<=n; v++) {
			fprintf (stderr, "%4d:", v);
			SPGArc *a, *end;
			for (GetBounds(v,a,end); a<end; a++) {
				fprintf (stderr, " %d:%d:%d ", a->head, a->cost, a->label);
			}
			fprintf (stderr, "\n");
		}		
	}


	void ReadSTP (char *file) {
		const int BUFSIZE = 2048;
		char buffer [BUFSIZE];
		int v, w;
		double cost;
		int n, m, t, tcount, ecount;
		double x, y;
		t = n = m = -1;
		tcount = ecount = 0;
		int tdeclared = -1;
		double mincost = 0.0000001;
		double fixed; //total fixed cost
		bool FOUND_PRESOLVE=true;
		
        ifstream FIC;
        FIC.open(file);
        if (FIC.fail())
        {
            printf("Fail to open the input file %s\n", file);
            return;
        }else{
			printf("read instances begin: %s!\n",file);
		}
        
		while (FIC.getline(buffer,BUFSIZE)) 
        {			
			//cout<<buffer<<endl;
			if (FOUND_PRESOLVE) {
				if (sscanf(buffer, "Fixed %lg\n", &fixed) == 1) {
					fprintf (stderr, "Read fixed cost: %.3f\n", fixed);
					break;////ignore fixed
					gd.IncFixed(fixed);
				}
				//continue; //ignore edges and terminals in presolve section
			} 

			if ((sscanf (buffer, "E %d %d %lf", &v, &w, &cost)>0) || (sscanf (buffer, "A %d %d %lf", &v, &w, &cost)>0)) {
				ecount ++; 
				if (cost == 0) {
					cost = mincost;
				}
				AddEdge(v,w,(EdgeCost)cost);
			} else if (sscanf (buffer, "T %d", &v)>0) {
				tcount ++; 
				MakeTerminal(v);
			} else if (sscanf (buffer, "Terminals %d", &t)>0) {
				tdeclared = t;
			} else if (sscanf (buffer, "Nodes %d", &n)>0) {
				SetVertices(n);
			} else if (sscanf (buffer, "Edges %d", &m)>0) {
				SetEdges(m);
			}
		}
		FIC.close();

		fprintf (stderr, "Done reading compact description of instance (n=%d m=%d t=%d).\n", gd.n, gd.GetEdgeCount(), gd.tcount); fflush(stderr);
		Commit();
	}

	inline EdgeCost GetCost (int e) const {
		return gd.edges[e].cost;
	}
    inline EdgeCost GetMinCost() {
		EdgeCost mincost = GetCost(1);
		for (int e=2; e<=m; e++) {
			if (GetCost(e) < mincost) mincost = GetCost(e);
		}
		return mincost;
	}

	inline int GetFirstEndpoint(int e) const {return gd.edges[e].v;}
	inline int GetSecondEndpoint(int e) const {return gd.edges[e].w;} 
	inline void GetEndpoints (int e, int &v, int &w) {
		v = gd.edges[e].v;
		w = gd.edges[e].w;
	}

	inline void GetDirectedEndpoints(int a, int &v, int &w) {
		if (a > m) {
			a -= m;
            v = gd.edges[a].w;
            w = gd.edges[a].v;
        } else {
			v = gd.edges[a].v;
            w = gd.edges[a].w;
        }
    }

	inline int GetOther (int e, int v) {
		int a, b;
		GetEndpoints(e,a,b);
		return (v==a) ? b : a;
	}



    /*for partition*/
	void ReadSubSTP (Graph &gg,int flagbelong) {

		int v, w;
		double cost;
		int n, m,nn=0;
		double x, y;
		n = m = 0;
		int n0=0;

		for(int v=1;v<=gg.VertexCount();v++)
		if(gg.GetBelong_to_Patition(v)==flagbelong)
		{
			nn++;
		}
		SetVertices(nn);
		//SetVertices(gg.VertexCount());
	
		for(int v=1;v<=gg.VertexCount();v++)
		if(gg.GetBelong_to_Patition(v)==flagbelong)
		{
			n++;
			newindex[v]=n;
			oldindex[n]=v;
			if(gg.IsTerminal(v))
			{
				MakeTerminal(n);
				//MakeTerminal(v);
			}
		}

		for(int i=1;i<=gg.EdgeCount();i++)
		{
			gg.GetArcEndpoints(i,v,w);
			if(gg.GetBelong_to_Patition(v)!=flagbelong||gg.GetBelong_to_Patition(w)!=flagbelong)
				continue;
			cost=gg.GetCost(i);
			AddEdge(newindex[v],newindex[w],(EdgeCost)cost);
			//AddEdge(v,w,(EdgeCost)cost);
            m++;
			oldedgeid[m]=i;
		}
        SetEdges(m);

		Commit();
	}


	void ReadSubSTP (Graph &gg,int flagbelong,vector<int>zeroedge) {

		int v, w;
		double cost;
		int n, m,nn=0;
		double x, y;
		n = m = 0;
		int n0=0;

		for(int v=1;v<=gg.VertexCount();v++)
		if(gg.GetBelong_to_Patition(v)==flagbelong)
		{
			nn++;
		}
		SetVertices(nn);
		//SetVertices(gg.VertexCount());
	
		for(int v=1;v<=gg.VertexCount();v++)
		if(gg.GetBelong_to_Patition(v)==flagbelong)
		{
			n++;
			newindex[v]=n;
			oldindex[n]=v;
			if(gg.IsTerminal(v))
			{
				MakeTerminal(n);
				//MakeTerminal(v);
			}
		}

		for(int i=1;i<=gg.EdgeCount();i++)
		{
			gg.GetArcEndpoints(i,v,w);
			if(gg.GetBelong_to_Patition(v)!=flagbelong||gg.GetBelong_to_Patition(w)!=flagbelong)
				continue;
			cost=gg.GetCost(i);
			if(zeroedge[cost]==1)
				cost=0.00000001;
				
			AddEdge(newindex[v],newindex[w],(EdgeCost)cost);
			//AddEdge(v,w,(EdgeCost)cost);
            m++;
			oldedgeid[m]=i;
		}
        SetEdges(m);

		Commit();
	}


	void ReadSubSTP (Graph &gg,vector<EdgeDescriptor>edgetoadd) {

		int v, w;
		double cost;
		int n, m;
		double x, y;
		n = m = 0;

		SetVertices(gg.VertexCount());
	
		for(int v=1;v<=gg.VertexCount();v++)
		{
			if(gg.IsTerminal(v))
			{
				MakeTerminal(v);
			}
		}

		for(int i=1;i<=gg.EdgeCount();i++)
		{
			gg.GetArcEndpoints(i,v,w);
			cost=gg.GetCost(i);
			AddEdge(v,w,(EdgeCost)cost);
            m++;
		}
		for(int i=0;i<edgetoadd.size();i++)
		{
			v=edgetoadd[i].v;
			w=edgetoadd[i].w;
			cost=edgetoadd[i].cost;
			AddEdge(v,w,(EdgeCost)cost);
            m++;
		}
        SetEdges(m);

		Commit();
	}

	
	/*for partition*/
};
