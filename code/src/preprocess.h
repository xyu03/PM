#pragma once
#include<string>
#include<cstring>
#include<cmath>
#include<cstdio>
#include<algorithm>
using namespace std;
class Preprocessing {
	public:
	template <class GRAPH> static void OutputSTP (const string &filename, GRAPH &g, vector<bool> &keep) {
		FILE *file = fopen (filename.c_str(), "w");
		int m = g.EdgeCount();
		int n = g.VertexCount();
		int newm = 0;

		for (int e=1; e<=m; e++) {
			if (keep[e]) newm++;
		}
		
		fprintf (file, "33d32945 STP File, STP Format Version 1.00\n");
		fprintf (file, "Section Comment\n");
		fprintf (file, "End\n\n");

		fprintf (file, "Section Graph\n");
		fprintf (file, "Nodes %d\n", n);
		fprintf (file, "Edges %d\n", newm);
		for (int e=1; e<=m; e++) {
			if (!keep[e]) continue;
			int v,w;
			g.GetEndpoints(e,v,w);
			EdgeCost cost = g.GetCost(e);
			fprintf (file, "E %d %d %.0f\n", v, w, (double)cost);
		}
		fprintf (file, "End\n\n");

		fprintf (file, "Section Terminals\n");
		fprintf (file, "Terminals %d\n", g.TerminalCount());
		for (int v=1; v<=n; v++) {
			if (g.IsTerminal(v)) fprintf (file, "T %d\n", v);
		}
		fprintf (file, "End\n\n");

		fclose(file);
	}

	template <class GRAPH> static void PruneNeighbors(GRAPH &g, int v, vector<bool> &arcs, vector<int> &degree, RFWStack<int> &candidates) {
		SPGArc *a, *end;
		for (g.GetBounds(v,a,end); a<end; a++) {
			if (!arcs[a->label]) continue;
			arcs[a->label] = false;
			int w = a->head;
			if (--degree[w] <= 2) {
				if (!g.IsTerminal(w)) candidates.push(w);
			}
		}
		degree[v] = 0;
	}

	static void DynamicToStatic (DynamicGraph &oldg, Graph &newg) {
		int newn = 0;
		int newm = 0;
		int oldn = oldg.VertexCount();
		const bool verbose = false;

		vector<int> old2new(oldn+1, -1); 
		if (verbose) fprintf (stderr, "Creating static graph (from %d)...\n", oldn); fflush(stderr);
		for (int v=1; v<=oldn; v++) {
			if (!oldg.IsTerminal(v) && (oldg.GetDegree(v) == 0)) continue;
			old2new[v] = ++newn;
		}
		if (verbose) fprintf (stderr, "New graph should have %d vertices.\n", newn);

		newg.SetVertices(newn);
		for (int v=1; v<=oldn; v++) {
			if (old2new[v]<0) continue;
			double x, y;
			if (oldg.IsTerminal(v)) newg.MakeTerminal(old2new[v]);
			SPGArc *a, *end;
			for (oldg.GetBounds(v,a,end); a<end; a++) {
				int w = a->head;
				if (w <= v) continue;
				newg.AddEdge(old2new[v],old2new[w],a->cost);
			}
		}

		newg.SetFixedCost(oldg.GetFixedCost());
		newg.Commit();
	
	}


        /// <summary> 
        /// Find MST of the distance network associated with a given set of bases.
        /// We can use uf to force the algorithm to think of two bases as a single one;
        /// uf will be altered during the algorithm, as different regions are joined.
        /// (Note: the algorithm could be implemented without uf---by coloring---but
        /// it's simpler to just do everything with uf.)
        /// </summary>
        /// <param name="solution">Output: the solution.</param>
        /// <param name="voronoi">Current voronoi diagram.</param>
        /// <param name="uf">List of regions; a voronoi region must be entirely within the same uf region.</param>
        /// <param name="pertcost">Perturbed edge costs (may be null).</param>

	static void BoruvkaGraph(Graph &g, SteinerSolution &solution, VoronoiData &voronoi, UnionFind &uf, EdgeCost *pertcost) {
		int v, n = g.VertexCount();
		int m = g.EdgeCount();
		EdgeCost solvalue = 0; //needed
        const bool verbose = false;

		//remember vertices with no outgoing boundary edges
		bool *hasneighbors = new bool [n+1];

		//count boundary regions in the current diagram
        int nregions = 0;
        for (v=1; v<=n; v++) {
			hasneighbors[v] = true; //as far as we know, everybody has a neighbor
            if (uf.Find(voronoi.GetBase(v)) == v) nregions ++;
		}

		//Boruvka's algorithm
        int *minarc = new int[n+1]; //minimum outgoing edge from the region based in v
        EdgeCost *minvalue = new EdgeCost[n+1]; //value of this neighbor

		int fakenumber = 0;
        int rounds = 0;
        bool changes = true;

		bool *useless = new bool [m+1];
		for (int e=1; e<=m; e++) {useless[e] = false;}

		const bool debug = false;
		while (changes && nregions > 1) {
			rounds ++;
			//if (rounds > 0) break;
			if (debug && rounds>4) break;
			if (verbose) fprintf (stderr, "%d:%d ", rounds, nregions);
			changes = false; //will be true if there is any progress in this round

			//initially, we don't know what are the arcs out of each component
			for (v=1; v<=n; v++) {
				minarc[v] = -1;
				minvalue[v] = 0;
			}

			int skipped = 0;

			//look for cheapest edge out of each component
			for (v=1; v<=n; v++) {
				if (!hasneighbors[v]) {
					skipped ++;
					if (verbose) fprintf (stderr, ".");
					continue; //skip internal vertices
				} else {
					if (verbose) fprintf (stderr, "|");
				}
				int bv = uf.Find(voronoi.GetBase(v)); //component containing v
				EdgeCost vdist = voronoi.GetDistance(v); //distance from v to its base
				int neighcount = 0;

				hasneighbors[v] = false; //unless we find an actual edge, this is false
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int w = a->head;					
					if (w>v) continue; //look at each edge once (COULD SKIP IF GRAPH GUARANTEED TO BE SORTED)
					if (useless[a->label]) {continue;}

					int bw = uf.Find(voronoi.GetBase(w)); //component containing w
					if (bv==bw) {
						useless[a->label] = true;
						continue; //no longer a boundary edge!
					}

					if (debug) {
						fakenumber += bw;
						changes = true;
						hasneighbors[v] = true;
						continue;
					}


					neighcount ++;
					hasneighbors[v] = true; //found a neighbor! (ACTUALLY, THIS SHOULD ONLY BE TRUE IF WE FIND AT LEAST 2, no?)
					//changes = true;
					//continue;

					//compute cost of corresponding boundary path
					int label = a->label;
					EdgeCost cost = (pertcost!=NULL) ? pertcost[label] : a->cost;
                    cost += vdist + voronoi.GetDistance(w);

					//update if better for bv
					if (minarc[bv]==-1 || (cost<minvalue[bv])) {minarc[bv]=label; minvalue[bv]=cost;}

					//update if better for bw
					if (minarc[bw]==-1 || (cost<minvalue[bw])) {minarc[bw]=label; minvalue[bw]=cost;}
				}

				// WARNING: THIS LOOKS VERY WRONG; I GUESS BOTH VERTICES COULD HAVE A SINGLE NEIGHBOR, AND NOBODY WILL BE VISITED. 
				// SOME TIEBREAKING IS NEEDED
				const bool AGGRESSIVE_OPTIMIZATION = false; //true;
				if (AGGRESSIVE_OPTIMIZATION) {
					if (neighcount <= 1) hasneighbors[v] = false; //if this vertex has a single neighbor, it will be merged with it in this iteration
				}
			}

			if (debug && fakenumber==314123411) {fprintf (stderr, "Miracle!\n");}

			//join each region to its best neighbor
			int bvisited = 0;
			int b;
			for (b=1; b<=n; b++) {
				if (minarc[b] >= 0) { //best neighbor defined...
					bvisited ++;
					int w = 0;
					g.GetEndpoints(minarc[b], v, w);
					int bv = voronoi.GetBase(v);
					int bw = voronoi.GetBase(w);
					if (uf.Find(bv) != uf.Find(bw)) { //join regions if not already joined [COULD UNITE DIRECTLY AND CHECK BOOL]
						changes = true;
						uf.Union(bv, bw);
						nregions--;
						solution.Insert(minarc[b]);

						// WARNING: THIS IS NOT STOPPING AS SOON AS IT SHOULD
						while (v != bv) {
							int f = voronoi.GetParentArc(v);
							if (!solution.Insert(f)) break;//; // fprintf (stderr, ".");
							v = g.GetOther(f, v);
						}
						while (w != bw) {
							int f = voronoi.GetParentArc(w);
							if (!solution.Insert(f)) break;
							w = g.GetOther(f, w);
						}
					}
				}
			}
			//fprintf (stderr, "s%.3f b%.3f r%d   ", (double)skipped / (double)n, (double)bvisited /(double)n, nregions);
		}
		//fprintf (stderr, "\n");


		if (verbose) fprintf(stderr, "\n"); //Done computing shortest paths; solution value so far is %d.\n", solvalue);
            //Console.Error.WriteLine("rounds:{0}", rounds);
		delete [] minvalue;
		delete [] minarc;
		delete [] hasneighbors;
		delete [] useless;
	}


	static void CopyFromStaticGraph (DynamicGraph &dg, Graph &sg, vector<bool> &arcs) {
		const bool verbose = false;
		if (verbose) fprintf (stderr, "Copying from static graph.\n");

		dg.SetVertices(sg.VertexCount());
		dg.SetEdges(sg.EdgeCount());
		for (int e=1; e<=sg.EdgeCount(); e++) {
			int v, w;
			sg.GetEndpoints(e,v,w);
			dg.BatchAddEdge(v,w,sg.GetCost(e));
		}
		for (int v=1; v<=sg.VertexCount(); v++) {
			if (sg.IsTerminal(v)) dg.MakeTerminal(v);
		}
		dg.SetFixedCost(sg.GetFixedCost());
		dg.Commit();

		// WARNING: THIS SHOULD BE BATCHED
		int hidden = 0;
		if (verbose) fprintf (stderr, "Hiding edges that have already been deleted!\n");
		for (int e=1; e<=sg.EdgeCount(); e++) {
			if (!arcs[e]) dg.HideEdge(e);
		}
		if (verbose) fprintf (stderr, "done hiding.\n");
		if (hidden>0) fprintf (stderr, "h%d ", hidden);
	}

	static void OneTwoPrep (Graph &sg, vector<bool> &arcs) {
		DynamicGraph g;
		CopyFromStaticGraph (g, sg, arcs);
		int m = g.EdgeCount();
		int n = g.VertexCount();
		int v, w;

		const int VERBOSE_LEVEL = 1;
		bool DO_NOTHING = false;
		const bool BOTTLE_TEST = true;
		const bool verbose = false;

		int MAX_DEGREE = 10;


		if (VERBOSE_LEVEL>=1) fprintf (stderr, "OneTwoPrep(%d,%d) ", n, m);
		//fprintf (stderr, "Vertex %d has degree %d\n", DEBUG_VERTEX, g.GetDegree(DEBUG_VERTEX));
		EdgeCost infinity = 10e32;
		vector<EdgeCost> distance (n+1, infinity);

		//vector<int> edgelist; // stack of edges to process
		vector<int> onelist;
		vector<int> twolist;
		vector<int> bottlelist; 
		// first copy all edges incident vertices of degree 1 and 2
		for (int e=1; e<=m; e++) {
			if (!arcs[e]) continue;
			g.GetEndpoints(e,v,w);
			if (g.GetDegree(v)==1 || g.GetDegree(w)==1) {onelist.push_back(e); continue;}
			if (g.GetDegree(v)<=2 || g.GetDegree(w)<=2) {twolist.push_back(e); continue;}
			if (BOTTLE_TEST) bottlelist.push_back(e);
		}

		while ((onelist.size()>0 || twolist.size()>0 || bottlelist.size()>0) && !DO_NOTHING) {
			int e;
			if (onelist.size() > 0) {
				e = onelist.back();
				onelist.pop_back();
			} else if (twolist.size() > 0) {
				e = twolist.back();
				twolist.pop_back();
			} else if (bottlelist.size() > 0) {
				//fprintf (stderr, "B");
				e = bottlelist.back();
				bottlelist.pop_back();
			}
			if (!arcs[e]) continue;

			g.GetEndpoints(e,v,w);
			if (v==w) {
				if (VERBOSE_LEVEL>=3) fprintf (stderr, "FOUND POTENTIAL SELF-LOOP!");
				if (g.IsTerminal(v)) {
					fprintf (stderr, "Apparently %d is a selfloop.\n", v);
				}
				arcs[e] = false;
				g.HideEdge(e);
				continue;
			}

			if ((g.IsTerminal(v) && g.GetDegree(v)>1) || g.GetDegree(v) > g.GetDegree(w)) std::swap(v,w);

			// v has the lowest degree / is the non-terminal


		
			// v is the vertex we will eliminate
			if (g.GetDegree(v)==1) {
				if (verbose) fprintf (stderr, "Removing (%d,%d) with degrees (%d,%d) and length %.0f.\n", v, w, g.GetDegree(v), g.GetDegree(w), g.GetCost(e));
				arcs[e] = false;
				g.HideEdge(e);

				
				if (g.IsTerminal(v)) {
					g.IncFixedCost(g.GetCost(e));
					g.UnmakeTerminal(v);
					if (!g.IsTerminal(w)) g.MakeTerminal(w);
					if (VERBOSE_LEVEL>=2) fprintf (stderr, "Removed %d, made %d terminal (%d).\n", v, w, g.GetDegree(w));
				}


				if (g.GetDegree(w)==2) {
					twolist.push_back(g.GetFirstEdgeLabel(w));
				} else if (g.GetDegree(w)==1) {
					onelist.push_back(g.GetFirstEdgeLabel(w));
				}

				if (g.GetDegree(w)==0 && g.IsTerminal(w)) {
					int vleft = 0;
					for (int y=1; y<=g.VertexCount(); y++) {
						if (g.GetDegree(y)>0) vleft ++;
					}
					fprintf (stderr, "neighbor (%d) of degree-one vertex (%d) was a degree-one terminal; vertices left=%d.\n", w, v, vleft);
				}

				// should already be in bottlelist, I guess
				continue;
			}

			if (g.IsTerminal(v)) continue;

			if (g.GetDegree(v)==2) {
				int neighbors[2];
				int edges[2];
				int ncount = 0;
				EdgeCost length = 0;

				// get list of neighbors
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					if (ncount == 2) {
						fprintf (stderr, "Vertex %d already has %d %d as neighbors.\n", v, neighbors[0], neighbors[1]);
					}
					neighbors[ncount] = a->head;
					edges[ncount] = a->label;
					length += a->cost;
					arcs[a->label] = false;
					ncount ++;
				}


				g.HideEdge(edges[0]);
				g.HideEdge(edges[1]);


				if (neighbors[0] == neighbors[1]) {
					if (verbose) fprintf (stderr, "Removing potential self-loop!\n");
					int u = neighbors[0];
					if (g.GetDegree(u)==1) onelist.push_back(g.GetFirstEdgeLabel(u));
					else if (g.GetDegree(u)==2) twolist.push_back(g.GetFirstEdgeLabel(u));
					else if (BOTTLE_TEST) bottlelist.push_back(g.GetFirstEdgeLabel(u));
				} else {
					if (verbose) fprintf (stderr, "Removed <%d-%d-%d> with length %.8f + %.8f.\n", neighbors[0], v, neighbors[1], g.GetCost(edges[0]), g.GetCost(edges[1]));

					int eid = g.InsertEdge(neighbors[0], neighbors[1], length);

					if (eid != (int)arcs.size()) {
						fprintf (stderr, "New id is %d, current size is %d.\n", eid, arcs.size());
					}
					arcs.push_back(true);

					int x = neighbors[0];
					int y = neighbors[1];

					//fprintf (stderr, "
					if (g.GetDegree(x)==1 && g.GetDegree(y)==1) onelist.push_back(eid);
					else if (g.GetDegree(x)==2 && g.GetDegree(y)==2) twolist.push_back(eid);
					else if (BOTTLE_TEST) bottlelist.push_back(eid);
				}
				continue;

			}

			// OK, we'll try to remove the edge by bottleneck distance
			if (!(g.IsTerminal(w) && g.GetDegree(w)==1)) {
				//fprintf (stderr, "B");

				if (g.GetDegree(v)+g.GetDegree(w) >= 2*MAX_DEGREE) {
					if (VERBOSE_LEVEL>=2) fprintf (stderr, "!");
					continue;

				}

				EdgeCost length = g.GetCost(e);

				SPGArc *a, *end;
				bool prune = false;

				// find distances from v
				for (g.GetBounds(v,a,end); a<end; a++) {
					if (arcs[a->label] && a->label!=e) {distance[a->head] = a->cost;}
				}

				if (distance[w] <= length) {
					if (VERBOSE_LEVEL>=2) fprintf (stderr, "Prunning at equality!\n");
					prune = true;
				} else {
					// try to find matches
					for (g.GetBounds(w,a,end); a<end; a++) {
						if (arcs[a->label] && a->label!=e) {
							if (distance[a->head] + a->cost <= length) {
								prune = true; break;
							}
						}
					}
				}

				// reset distances
				for (g.GetBounds(v,a,end); a<end; a++) {
					if (arcs[a->label]) {
						if (a->head == w) {
							if (a->cost > length) {
								if (BOTTLE_TEST) bottlelist.push_back(a->label); //CAREFUL: MAY ADD SEVERAL TIMES
							}
						}
						distance[a->head] = infinity;
					}
				}

				if (prune) {
					if (VERBOSE_LEVEL>=2) fprintf (stderr, "REMOVING EDGE (%d,%d)\n", v, w);
					arcs[e] = false;
					g.HideEdge(e);
					if (g.GetDegree(v) == 2) twolist.push_back(g.GetFirstEdgeLabel(v));
					else if (g.GetDegree(v) == 1) onelist.push_back(g.GetFirstEdgeLabel(v));
					if (g.GetDegree(w) == 2) twolist.push_back(g.GetFirstEdgeLabel(w));
					else if (g.GetDegree(w) == 1) onelist.push_back(g.GetFirstEdgeLabel(w));
				}

			}

			//fprintf (stderr, "Should not have gotten here.");
		}
		
		Graph newg;
		DynamicToStatic(g, newg);
		sg = newg;
	}

	// recursively remove corners
	static void PrepCorners(Graph &sg, vector<bool> &arcs, vector<int> &degree) {
		DynamicGraph g;
		CopyFromStaticGraph (g, sg, arcs);
		
		int m = g.EdgeCount();
		int n = g.VertexCount();
		int v, w;



		fprintf (stderr, "Preprocessing corners (%d,%d)!\n", n, m);

		EdgeCost infinity = 10e32;
		vector<EdgeCost> distance (n+1, infinity);


		//reset degrees
		for (int v=1; v<=n; v++) {
			degree[v] = 0;
		}

		//compute actual degrees
		for (int e=1; e<=m; e++) {
			if (arcs[e]) {
				g.GetEndpoints(e,v,w);
				fprintf (stderr, "(%d,%d)", v, w);
				degree[v]++;
				degree[w]++;
			}
		}

		RFWStack<int> candidates(n);
		for (int v=1; v<=n; v++) {
			if (degree[v] <= 2 && !g.IsTerminal(v)) candidates.push(v);
		}

		while (!candidates.isEmpty()) {
			v = candidates.pop();
			fprintf (stderr, "Checking vertex %d with degree %d.\n", v, degree[v]);

			if (degree[v] == 1) {
				PruneNeighbors(g, v, arcs, degree, candidates);
				continue;
			}

			if (degree[v] == 2) {
				fprintf (stderr, "Vertex %d has degree %d.\n", v, degree[v]);
				int neighbor[2];
				EdgeCost cost[2];
				int ncount =0;
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					if (arcs[a->label]) {
						cost[ncount] = a->cost;
						neighbor[ncount] = a->head;
						if (++ncount == 2) break;
					}
				}
				fprintf (stderr, "Vertex %d has neighbors %d and %d.\n", v, neighbor[0], neighbor[1]);
				EdgeCost length = cost[0] + cost[1];
				bool prune = false;
				if (neighbor[0] == neighbor[1]) {
					prune = true;
				} else {
					for (g.GetBounds(neighbor[0],a,end); a<end; a++) {
						if (arcs[a->label] && a->head!=v) {distance[a->head] = a->cost;}
					}

					for (g.GetBounds(neighbor[1],a,end); a<end; a++) {
						if (arcs[a->label] && a->head!=v) {if (distance[a->head] + a->cost <= length) prune = true;}
					}

					for (g.GetBounds(neighbor[0],a,end); a<end; a++) {
						if (arcs[a->label] && a->head!=v) {distance[a->head] = infinity;}
					}
				}

				if (prune) {
					PruneNeighbors(g, v, arcs, degree, candidates);
				}
			}
		}
	}

	//-----------------------------------------------
	// 移除度为1的点
	//-----------------------------------------------
	static void PrepDegreeOne(Graph &g, vector<bool> &arcs, vector<int> &degree) {
		int n = g.VertexCount();

		// warning: could loop over edges instead
		int v;
		for (v = 1; v <= n; v++) {
			degree[v] = 0;
			SPGArc *a, *end;
			for (g.GetBounds(v,a,end); a<end; a++) {
				if (arcs[a->label]) degree[v]++;
			}
		}

		int removed = 0;
		int terminals = 0;
		int i;
		for (i=1; i<=n; i++) {
			v = i;
			while (!g.IsTerminal(v) && degree[v]==1) {
				SPGArc *a, *end;
				for (g.GetBounds(v,a,end); a<end; a++) {
					int e = a->label;
					if (arcs[e]) {
						removed ++;
						arcs[e] = false;
						degree[v] --;
						int w = a->head;
						degree[w] --;
                        if (degree[w] == 1) {v = w;}
						break;
					}
				}
			}
            
			if (degree[v] == 1) {
				terminals++;
			}
		}
		fprintf (stderr, "Removed %d vertices, and there are still %d terminals of degree 1.", removed, terminals);
	}


	/*-------------------------------------------------------------
	 | 创建一个仅包含被保留边的新的图
	 *------------------------------------------------------------*/
	static void CommitSubgraph (Graph &newg, Graph &g, vector<bool> &edgekeep) {
		int oldn = g.VertexCount();
		int oldm = g.EdgeCount();
		int newm = 0;
		int v, w, e;
		const bool verbose = false;

		if (verbose) fprintf (stderr, "Creating subgraph after preprocessing.\n");

		vector<int> degree(oldn+1, 0);
		
		// figure out the degree
		for (e=1; e<=oldm; e++) {
			if (edgekeep[e]) {
				g.GetEndpoints(e, v, w);
				degree[v] ++;
				degree[w] ++;
				newm ++;
			}
		}

		int newcount = 0;
		vector<int> old2new(oldn+1, -1); 
		int onecount = 0;
		int twocount = 0;
		for (v=1; v<=oldn; v++) {
			if (degree[v] == 0) continue;
			if (degree[v] == 1 && !g.IsTerminal(v)) {onecount ++;} // continue;}
			if (degree[v] == 2 && !g.IsTerminal(v)) {twocount ++;}
			old2new[v] = ++newcount;
		}

		if (verbose) fprintf (stderr, "Skipped %d vertices of degree 1 and counted %d of degree 2; graph has %d vertices.\n", onecount, twocount, newcount);
		//fprintf (stderr, "Graph has %d vertices.\n", newcount);

		// FINALLY CREATE NEW GRAPH
		newg.SetVertices(newcount);
		newg.SetEdges(newm);
		newg.SetFixedCost(g.GetFixedCost());
		for (e=1; e<=oldm; e++) {
			if (edgekeep[e]) {
				g.GetEndpoints(e, v, w);
				newg.AddEdge(old2new[v],old2new[w],g.GetCost(e));
			}
		}
		for (v=1; v<=oldn; v++) {
			if (old2new[v] < 0) continue;
			if (g.IsTerminal(v)) newg.MakeTerminal(old2new[v]);
		}
		newg.Commit();

		if (verbose) fprintf (stderr, "Created actual graph with %d vertices and %d edges.\n", newg.VertexCount(), newg.EdgeCount());

	}

	static void RunPreprocessing (Graph &sg, const bool runOneTwo) {
		
		int origm = sg.EdgeCount();
		double threshold = 0.01;
		int rounds;
		for (rounds=0; rounds<999; rounds++) {
			int oldedges = sg.EdgeCount();
			PrepBottleneck(sg);
			int m = sg.EdgeCount();		

			if (runOneTwo) {
				vector<bool> keep(m + 1, true);
				OneTwoPrep(sg, keep);
			}
			int newedges = sg.EdgeCount();
			if (newedges == 0) {
				fprintf (stderr, "Instance solved to optimality!\n");
				break;
			}
			double reduction = (double)(oldedges - newedges) / (double)oldedges; 
			fprintf (stderr, "R%d:%.3f ", rounds, reduction);
			if (reduction < threshold) {
				fprintf (stderr, "!\n"); //"Reduction was only %.3f; stopping.\n", reduction);				
				break;
			}
		}
		
		printf ("Done preprocessing.\n");
		printf ("preprounds %d\n", rounds+1);
		printf ( "prepthreshold %.3f\n", threshold);
		
		double reduction = (origm - sg.EdgeCount()) / (double)origm;
		printf ( "prepreduction %.6f\n", reduction); 
		printf ("prepreductionpct %.6f\n", 100.0 * reduction); 
		printf ("prepsolved %d\n", sg.EdgeCount()==0);
		printf ("prepvertices %d\n", sg.VertexCount());
		printf ("prepedges %d\n", sg.EdgeCount());
		printf ("prepterminals %d\n", sg.TerminalCount());
		printf ("prepfixedcost %.6f\n",sg.GetFixedCost());

	}

	//**瓶颈距离：（u,v），如果存在一条u-v瓶颈路径B(u,v)>=e(u,v),删除(u,v)
	//计算出基于terminal的Voronoi diagram图，然后利用Boruvka最小生成树算法将这些块连接起来
	static void PrepBottleneck(Graph &g) {
		int n = g.VertexCount();
		int m = g.EdgeCount();
        int v;
        fprintf (stderr, "Preprocessing instance... ");

		const bool verbose = false;
		SteinerSolution solution(&g);

		UniverseSet baselist(n);

		for (v=1; v<=n; v++) {
			if (g.IsTerminal(v)) baselist.Insert(v);
			//degree[v] = g.GetDegree(v);
		}

		BinaryHeap <EdgeCost> binheap(n); // = new BinaryHeap<ArcCost> (n);
		VoronoiData vordata(n);// = new VoronoiData(n);
		UnionFind uf(n); // = new UnionFind(n);
		Basics::ComputeVoronoi(g,vordata,baselist,binheap,NULL); 
		BoruvkaGraph(g,solution,vordata,uf,NULL);

		if (verbose) fprintf (stderr, "Found solution worth %.3f.\n", solution.GetCost());

		// find maximum bottleneck length within the graph
		EdgeCost maxbottle = 0;
		for (int e=1; e<=m; e++) {
			if (!solution.Contains(e)) continue;
			int v, w;
			g.GetEndpoints(e,v,w);
			if (vordata.GetBase(v)!=vordata.GetBase(w)) { //found bridge
				EdgeCost pathcost = vordata.GetDistance(v) + vordata.GetDistance(w) + g.GetCost(e);
				if (pathcost > maxbottle) maxbottle = pathcost;
			}
		}
		
		// remove stupid edges
		int neliminated = 0;
		vector<int> elist;
		vector<bool> keep (m+1,true);
		vector<EdgeCost> priority (m+1,-1); 
		for (int e=1; e<=m; e++) {
			if (solution.Contains(e)) {
				int v,w;
				g.GetEndpoints(e,v,w);
				int bv = vordata.GetBase(v);
				int bw = vordata.GetBase(w);
				if (bv!=bw) {
					priority[e] = g.GetCost(e) + vordata.GetDistance(v) + vordata.GetDistance(w); 
				} 
			} else {
				EdgeCost cost = g.GetCost(e);
				priority[e] = cost; //min(min(dv,dw),cost); //to make sure it appears after all edges
				if (cost < maxbottle) continue;
				EdgeCost dv = vordata.GetDistance(g.GetFirstEndpoint(e));
				if (cost < dv) continue;
				EdgeCost dw = vordata.GetDistance(g.GetSecondEndpoint(e));
				if (cost < dw) continue;
				neliminated ++;
				keep[e] = false;
			}
		}
		for (int e=1; e<=m; e++) {
			if (keep[e] && priority[e]>=0) {
				elist.push_back(e);
			}
		}
		if (verbose) fprintf (stderr, "\n");
		int kept = (int)elist.size();

		//fprintf (stderr, "Maximum bottleneck length is %.3f.\n", maxbottle);
		if (verbose) fprintf (stderr, "Maxbotttle=%.3f eliminated %d/%d edges, keeping %d(%d).\n", maxbottle, neliminated, m, m - neliminated, kept);
		//fflush(stderr);
		//const bool verbose = false;

		bool FANCY_BOTTLENECK = true;
		if (FANCY_BOTTLENECK) {
			if (verbose) fprintf(stderr, "elist.size=%d, kept=%d, priority.size=%d, edge.size=%d\n", elist.size(), kept, priority.size(), solution.EdgeCapacity());
			
			// Sort all relevant edges by their bottleneck relevance. In case of ties, solution edges are inserted first.
			sort(elist.begin(), elist.begin() + kept, [&](int x, int y) {
				if (verbose) { 
					assert(x < priority.size()); 
					assert(y < priority.size()); 
					assert(x < solution.EdgeCapacity()); 
				} 
				return (priority[x] < priority[y] || (priority[x] == priority[y] && solution.Contains(x) && !solution.Contains(y)));
			});
			
			if (verbose) fprintf(stderr, "Sorted.\n");

			EdgeCost curbottle = 0;
			UnionFind newuf(n);
			///对于一条边(v,w)，如果两个端点不同时在解中则可以被删除
			///如果同时存在解T中，当删除这条边产生Tv和Tw两个子树，
			///如果Tv和Tw之间存在一个瓶颈距离小于e(v,w)，则这条边也可以被删除
			for (int i=0; i<kept; i++) {
				//fprintf (stderr, "-");
				int e = elist[i];
				int v, w;
				g.GetEndpoints(e,v,w);
				bool localdebug = false;

				if (solution.Contains(e)) {
					if (localdebug) fprintf (stderr, "Edge %d belongs to the solution.\n", e);
					int bv = vordata.GetBase(v);
					int bw = vordata.GetBase(w);
					newuf.Union(bv,bw);
					if (verbose && (curbottle!=priority[e])) fprintf (stderr, "   B%d ", priority[e]);
					curbottle = priority[e];
					//if (curbottle==2) exit(-1);
				} else {
					//this is an edge we have to check
					EdgeCost cost = priority[e];
					if (cost!=g.GetCost(e)) fprintf (stderr, ".");
					int v, w;
					g.GetEndpoints(e,v,w);
					if (vordata.GetParentArc(v)==e || vordata.GetParentArc(w)==e) {
						if (localdebug) fprintf (stderr, "Edge %d is a parent edge!\n", e);
						if (verbose) fprintf (stderr, ".%.0f", g.GetCost(e));
						keep[e] = true;
					} else {
						int bv = vordata.GetBase(v);
						int bw = vordata.GetBase(w);
						if (localdebug) fprintf (stderr, "Edge %d is a candidate!\n", e);
						if (newuf.Find(bv) == newuf.Find(bw)) {
							if (localdebug) fprintf (stderr, "Edge %d has both endpoints in the same component!\n", e);
							EdgeCost dv = vordata.GetDistance(v);
							EdgeCost dw = vordata.GetDistance(w);

							EdgeCost localbottle = max(dv,dw);
							if (cost >= localbottle) {
								
								if (DEBUG_VERTEX>=0) {
									if (v==DEBUG_VERTEX || w==DEBUG_VERTEX) {
										fprintf (stderr, "Removing edge (%d,%d).\n", v, w);
										fprintf (stderr, "dv=%.3f dw=%.3f cost=%.3f (v=%d w=%d bv=%d bw=%d) (e=%d pw=%d pv=%d) lb=%.3f cb=%.3f\n", dv, dw, cost, v, w, bv, bw, e, vordata.GetParentArc(v), vordata.GetParentArc(w), localbottle, curbottle);
									}
								}
								
								keep[e] = false;
							}
						}
					}
					//fprintf (stderr, "%d", (int)keep[e]);
				}
			}
		}

		int kcount = 0;
		for (int e=1; e<=m; e++) {
			if (keep[e]) {
				kcount++;
			}
		}
		if (verbose) fprintf (stderr, "\nKept %d edges eventually.\n", kcount); fflush(stderr);

		const bool REMOVE_NON_SP_EDGES = true;
		if (REMOVE_NON_SP_EDGES) {
			// remove edges with 
			int removed = 0;
			vector<int> last(n+1,-1);
			vector<EdgeCost> value(n+1,0);
			for (int e=1; e<=m; e++) {
				if (!keep[e]) continue;
				EdgeCost ecost = g.GetCost(e);
				int v,w;
				g.GetEndpoints(e,v,w);
				for (int s=0; s<=1; s++) {
					int x = (s==0) ? v : w; 
					SPGArc *a, *end;
					for (g.GetBounds(x,a,end); a<end; a++) {
						if (!keep[a->label]) continue;
						if (e==a->label) continue;
						int y = a->head;
						EdgeCost cost = a->cost;
						if (s==0) {
							if (last[y]!=x || value[y]>cost) {
								last[y] = x;
								value[y] = cost;
							}
						} else {
							if ((last[y]==v) && (cost+value[y] <= ecost)) {
								if (verbose) fprintf (stderr, "Removing %d-%d-%d\n", v, y, w);
								keep[e] = false;
							}
						}
					}
				}
				if (!keep[e]) removed++;
			}
			fprintf (stderr, "Removed %d additional edges, left with %d.\n", removed, kcount - removed);
		}

		fprintf (stderr, "Removed non-sp edges.\n");
		fflush(stderr);


		vector<int> degree(n+1,0);
		int maxdegree = 0;
		for (int e=1; e<=m; e++) {
			if (keep[e]) {
				int dv = ++degree[g.GetFirstEndpoint(e)];
				int dw = ++degree[g.GetSecondEndpoint(e)];
				maxdegree = max(maxdegree,max(dv,dw));
			}
		}

		const bool REPORT_DEGREES = false;
		if (REPORT_DEGREES) {
			vector<int> degcount(maxdegree+1,0);
			for (int v=1; v<=n; v++) {
				degcount[degree[v]] ++;
			}
			fprintf (stderr, "Degrees:\n");
			for (int d=0; d<=maxdegree; d++) {
				if (degcount[d]>0) fprintf (stderr, "%3d: %5d\n", d, degcount[d]);
			}
			fprintf (stderr, "done reporting degrees\n");
			fflush(stderr);
		}

		bool OUTPUT_FILES = false;

		if (OUTPUT_FILES) {

			vector<bool> vkeep(n+1, true);
			for (int v=1; v<=n; v++) {
				if (degree[v] == 0) vkeep[v] = false;
			}

			fprintf (stderr, "About to draw...\n"); fflush(stderr);
			fprintf (stderr, "Drawn.\n"); fflush(stderr);
			OutputSTP("test.stp", g, keep);
		}

		Graph newg;
		CommitSubgraph (newg, g, keep);
		g = newg;

	}
};
