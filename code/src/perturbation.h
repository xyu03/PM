#pragma once
#include <algorithm>
#include "config.h"
#include "graph.h"
#include "elite.h"
class Perturbator {
	int PERTURBATION_MODE;
	double base;
	double range;
	double threshold;
	SteinerConfig config;

	public:
		void SetParameters (int n, SteinerConfig &_config) {
			config = _config;
			PERTURBATION_MODE = config.PERTURBATION_MODE;
			base = config.PERT_FLOOR;
			range = config.PERT_RANGE;
			threshold = 0;
			if (PERTURBATION_MODE == 2) {
				threshold = sqrt((double)n) / (double)n;
			} else if (PERTURBATION_MODE == 4) {
				threshold = sqrt((double)n) / (double)(2*n);
			} else {
				threshold = (log(n)/log(2)) / (double)n;
			}
		}

		double GetPerturbation(double r) {
			//fprintf(stderr, "r is %.10f, threshold=%.10f\n", r, threshold);
			//cerr << r << endl;
			switch (PERTURBATION_MODE) {
				case 0: return base + r * range; break; //expected value is 1.5
				case 1: return base + r + 1/(10.0*r); break; //NOT TRUE EXPECTATION
				case 2: return (r >= threshold) ? 1 + range * r : range/r; break;
				case 3: return base + range * r*r; break;
				case 4: return (r >= threshold) ? 1 + range * r : range/r; break;
				case 5: return base + range * (1.0 - r*r); break;
				case 6: return base + 1.5 * r * r; break; //expected value is 1.5
				case 7: return (r >= threshold) ? 1 + range * r : max(r/threshold, 0.0000001); break;
				case 8: return 1.0 + range * r;
			}
			return 0;
		}

		void ResetRange (double r) {
			if (config.PERT_EXTRA > 0) {
				range = config.PERT_RANGE + config.PERT_EXTRA * r;
			}
			/*
			if (VARIABLE_RANGE) {
				if (3*r < 1) range = .25;
				else if (3*r > 2) range = 1;
				else range = 4;
			}*/
		}
	};



class PerturbationTools {
public:
	
	static void VertexPerturbation (Graph &g, vector<EdgeCost> &pertcost, SteinerConfig &config) {
		int n = g.VertexCount();
		int m = g.EdgeCount();

		Perturbator p;
		p.SetParameters(n, config);
		p.ResetRange(1.0*(rand()%1000)/1000);

		static bool first = true;
		if (first) {
			fprintf (stderr, "USING PERTURBATION MODE %d\n", config.PERTURBATION_MODE);
			first = false;
		}
		vector<double> vpert(n+1);

		for (int v=0; v<=n; v++) {
			//fprintf(stderr, "<<<<");
			vpert[v] = p.GetPerturbation(1.0*(rand()%1000)/1000);
		}

		// each edge is perturbed by the average of its endpoints
		for (int e = 1; e <= m; e++) {
			int v, w;
            g.GetEndpoints(e, v, w);
            double p = (vpert[v] + vpert[w]) / 2.0;
			EdgeCost cost =g.GetCost(e) * (EdgeCost)p;
            pertcost[e] = cost;
		}
	}



	static void OldVertexPerturbation (Graph &g, vector<EdgeCost> &pertcost) {
		const bool debug = false;
		//fprintf (stderr, "vp");

		int n = g.VertexCount();
        int m = g.EdgeCount();
        int divisor = 1;

        EdgeCost mine = g.GetMinCost();
		while (mine > 1000) {
			divisor *= 10;
            mine /= 10;
		}

		bool UNIFORM = false;

		int t = g.TerminalCount();

		vector<EdgeCost> vpert(n+1);
		for (int v=0; v<=n; v++) {
			if (UNIFORM) {
				// it was 100,200
				vpert[v] = rand()%100+100; //it's 100,120 for local stuff
				//vpert[v] = random.GetInteger(100, 200); 100 + 20 * g.GetDegree(v));
			} else {
				double r = 1.0*(rand()%1000)/1000;
				//r = (1 - r*r);
				//vpert[v] = 100 + 100 * (r + 1/r);
				vpert[v] = 100*r + 10/r;
				//vpert[v] = 100 + 20*r + 1/r;

				//fprintf (stderr, "%.1f ", vpert[v]);
			}
		}

		for (int e = 1; e <= m; e++) {
			int v, w;
            g.GetEndpoints(e, v, w);
            bool average = true;

            EdgeCost p;
			if (average) {
				p = (vpert[v] + vpert[w]) / 2;
			} else {
				/*
                    ArcCost minp = (vpert[v] < vpert[w]) ? vpert[v] : vpert[w];
                    ArcCost maxp = (vpert[v] > vpert[w]) ? vpert[v] : vpert[w];
                    p = (3*minp + maxp) / 4;
					*/
			}

			EdgeCost cost = (g.GetCost(e) / divisor) * p;
            pertcost[e] = cost;
		}
	}

	static void AdaptivePerturbation (Graph &g, vector<EdgeCost> &pertcost, SolutionPool &elite) {

		int m = g.EdgeCount();
		int solcount = elite.Count();

		//fprintf (stderr, "There should be %d elite solutions.\n", solcount);
		//fflush (stderr);

		for (int e=1; e<=m; e++) pertcost[e] = 0;

		// count the number of times each edge is used
		for (int i=1; i<=solcount; i++) {
			SteinerSolution *s = elite.GetReference(i);
			for (int e=1; e<=m; e++) {
				//fprintf (stderr, "%d ", e);
				//fflush (stderr);
				if (s->Contains(e)) pertcost[e] ++;
			}
		}

		//fprintf (stderr, "Done computing multiplicities.\n");

		double factor = (90.0 / (double)solcount);

		int maxc = 0;
		for (int e=1; e<=m; e++) {
			int c = (int)pertcost[e];
			if (c > maxc) maxc = c;
			int l=100;
			int r=110 + (int)ceil((double)c * factor);
			double rd=1.0*(rand()%(r-l+1)+l);
			pertcost[e] = g.GetCost(e) * rd;
		}

		fprintf (stderr, "%d ", maxc, solcount);
	}

    static void InitPerturbation(Graph &g, vector<EdgeCost> &pertcost, SteinerConfig &config) {
		//double p = config->PERT_VERTEX;
		bool USE_VERTEX_PERTURBATION = (1.0*(rand()%1000)/1000 < config.PERT_VERTEX);
		//fprintf (stderr, "USING VERTEX PERTURBATION? %d\n", USE_VERTEX_PERTURBATION);

		//bool USE_VERTEX_PERTURBATION = config->PERT_VERTEX; //random.GetInteger(1,2)==1; //
        if (USE_VERTEX_PERTURBATION) {
			VertexPerturbation(g, pertcost,  config);
			return;
        }
		Perturbator p;
		p.SetParameters(g.VertexCount(), config); 
		p.ResetRange(1.0*(rand()%1000)/1000);

		int m = g.EdgeCount();
		for (int e=1; e<=m; e++) {
			//pertcost[e] = (base + range * random.GetDouble()) * g.GetCost(e);
			pertcost[e] = p.GetPerturbation(1.0*(rand()%1000)/1000) * g.GetCost(e);
		}
	}

	static void  PartitionPerturbation(Graph &g,vector<EdgeCost> &pertcost,vector<int> &edgeappear,SteinerConfig &config)
	{
		int m=g.EdgeCount();
		Perturbator p;
		p.SetParameters(g.VertexCount(), config); 
		p.ResetRange(1.0*(rand()%1000)/1000);

		int flag=1;//rand()%2;

		for(int e=1;e<=m;e++)
		{
			pertcost[e]=g.GetCost(e)*p.GetPerturbation(1.0*(rand()%1000)/1000);
			double pe=1.0;//(1.0*edgeappear[e])/(1.0*edgeappear[0]);
			if(flag)
			{
				if(edgeappear[e])
					pertcost[e]=pertcost[e]*(1.0-0.75*pe*(1.0*(rand()%10000)/10000));
			}
			else
			{
				if(edgeappear[e])
					pertcost[e]=pertcost[e]*(1.0+0.75*pe*(1.0*(rand()%10000)/10000));
			}
		}
	}
};

