
#pragma once

#include <cstdio>
#include <cassert>
#include "graph.h"

class SteinerSolution {
public:
	Graph *g;
private:
	vector<bool> edge;  //incident vector of edges in the solution
	vector<int> degree; //current degrees of all solution vertices
	int nleaves;        //number of nonterminal leaves
	EdgeCost cost;      //total solution cost

	///增加一条边后更新端点
	inline void ProcessInsertion(int v) {
		int d = ++degree[v];
        if (d<=2 && !g->IsTerminal(v)) {
			if (d==1) nleaves ++;
            else nleaves --;
		}
	}

	///删除一条边后更新端点
	inline void ProcessRemoval (int v) {
		int d = --degree[v];
        if (d<2 && !g->IsTerminal(v)) {
			if (d==0) nleaves --;
            else nleaves ++;
		}
	}


public:
	inline int LeafCount() const {return nleaves;}
	inline bool Contains(int e) const {  return edge[e]; }
	inline size_t EdgeCapacity() const { return edge.size(); }

	void Reset() {
		fill_n(edge.begin(), g->EdgeCount()+1, false);
		fill_n(degree.begin(), g->VertexCount()+1, 0);
		nleaves = 0;
		cost=g->GetFixedCost();
	}

	inline bool Insert (int e) {
		if (Contains(e)) return false;
		edge[e] = true;
		cost += g->GetCost(e);
		ProcessInsertion(g->GetFirstEndpoint(e));
		ProcessInsertion(g->GetSecondEndpoint(e));
		return true;
	}

	int EdgeCount() {
		int ecount = 0;
		int m = g->EdgeCount();
		for (int e=1; e<=m; e++) {
			if (Contains(e)) ecount ++;
		}
		return ecount;
	}

	void Output(FILE *file) {
		fprintf (stderr, "Outputting solution..."); 
		int m = g->EdgeCount();
		fprintf (file, "m %d\n", EdgeCount());
		fprintf (file, "c %.0f\n", GetCost());
		for (int e = 1; e<=m; e++) {
			if (!Contains(e)) continue;
			int v, w;
			g->GetEndpoints(e, v, w);
			fprintf (file, "e %d %d %.0f\n", v, w, g->GetCost(e));
		}
		fprintf (stderr, "done.\n");
	}

	void Output (char *prefix) {
		char filename[2048];
		sprintf (filename, "%s.%09.0lf.sol", prefix, (double)GetCost());
		FILE *file = fopen (filename, "w");
		if (!file) {
			fprintf (stderr, "Could not open <%s> for output. Will not output solution.\n", filename);
			return;
		}
		Output(file);
		fclose(file);
		fprintf (stderr, "Solution written to <%s>.\n", filename);
		fflush(stderr);
	}




	/// Remove edge e from the current solution. Returns true iff
    /// the edge was actually in the solution before. 
    inline bool Remove(int e) {
		if (edge[e]==false) return false;
		edge[e] = false;
		cost -= g->GetCost(e);
		ProcessRemoval(g->GetFirstEndpoint(e));
		ProcessRemoval(g->GetSecondEndpoint(e));
		return true;
	}


	inline int GetDegree (int v) const {return degree[v];}
	inline EdgeCost GetCost () const {return cost;}

	bool IsBetter (SteinerSolution *s) {
		return (GetCost() < s->GetCost() - EDGE_COST_PRECISION);
	}

	void UpdateCost() {
		const bool verbose = false;
		if (verbose) fprintf (stderr, "Cost updated from %d to ", cost);
		cost=g->GetFixedCost();
		int m = g->EdgeCount();
		for (int e=1; e<=m; e++) {
			if (Contains(e)) cost += g->GetCost(e);
		}
		if (verbose) fprintf (stderr, "%d.\n", cost);
	}



	SteinerSolution (Graph *_g) {
		g = _g;
		edge.resize(g->EdgeCount()+1);
		degree.resize(g->VertexCount()+1);
		Reset();
	}

	SteinerSolution (SteinerSolution *s) {
		g = s->g;
		edge.resize(g->EdgeCount()+1);
		degree.resize(g->VertexCount()+1);
		Reset();
		CopyFrom(s);
	}

	void CopyFrom(SteinerSolution *s) {
		if (g != s->g) {
			g = s->g;
		} else {
			Reset();
		}
		int m = g->EdgeCount();
		for (int e=1; e<=m; e++) {
			if (s->Contains(e)) Insert(e);
		}

		//fprintf (stderr, "Created solution with %d arcs.\n", this->EdgeCount());

        //foreach (int e in ((SteinerSolution)s).ElementEnumerator()) {
        //        Insert(e);
        //}
            //Console.Error.WriteLine("copy not implemented");   
	}


	/// <summary>
	/// Compute degree of difference between this solution and ds, given 
	/// by (union - intersection) / intersection. In particular, the result
	/// is zero iff the solutions are identical, and 1 if they have no edge
    /// in common.
	double GetDifference(SteinerSolution *s) {
            //SteinerSolution s = (SteinerSolution)ds;
            //int common = 0;
            //int scount = 0;

		int ucount = 0; //edges.Count(); //number of edges in the union
        int icount = 0; //number of edges in the intersection

		int m = g->EdgeCount();
		for (int e=1; e<=m; e++) {
			int count = (Contains(e) ? 1 : 0) + (s->Contains(e) ? 1 : 0);
			if (count == 0) continue;
			//fprintf (stderr, "%d%d ", Contains(e), s->Contains(e));
			ucount ++; //1 or 2
			if (count == 2) {
				icount ++;
			}
		}

		//fprintf (stderr, "[%d/%d] ", icount, ucount);

        return (1.0 - (double)icount/(double)ucount);
	}

};

