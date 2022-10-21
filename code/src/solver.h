#include "elite.h"
#include "perturbation.h"

int MST (Graph &g, SteinerSolution &solution, UniverseSet &svertices) 
{
    bool verbose = false;
    int n = g.VertexCount();
    BinaryHeap<EdgeCost> heap(n);
    vector<int> parc (n+1);
    int r = Basics::PickRandomTerminal(g);
    parc[r] = 0;
    int nscanned = 0;

    //run Prim's algorithm
    solution.Reset();
    heap.Insert(r, 0);
    while (!heap.IsEmpty()) {
        unsigned int v;
        EdgeCost acost;
        heap.RemoveFirst(v,acost);
        if (v!=r) solution.Insert(parc[v]); //add edge (p(v),v) to solution

        nscanned ++;
        SPGArc *a, *end;
        for (g.GetBounds(v,a,end); a<end; a++) {
            int w = a->head; 
            if (!svertices.Contains(w)) continue; //we only care about svertex
            if (solution.GetDegree(w) > 0) continue; //vertex already in the new tree
            if (heap.Insert(w, a->cost)) parc[w] = a->label;
        }
    }
    return nscanned;
}


static bool MSTPrune(Graph &g, SteinerSolution &solution) {
		//return 0;
		EdgeCost original = solution.GetCost();
		int n = g.VertexCount();
		UniverseSet svertices(n); 
		Basics::MarkSolutionNodes(g, solution, svertices);
		MST(g,solution,svertices);
		Basics::Prune(g,solution);
		return (solution.GetCost() < original) ? 1 : 0;
	}

void RunLocalSearch (Graph &g, SteinerSolution &solution, int maxrounds, SteinerConfig &config, ExecutionLog &executionLogPtr) {
	bool RUN_Q = false;
	bool RUN_V = false;
	bool RUN_U = false;
	bool RUN_K = false;
	bool RESTRICT_K = false; //
	bool verbose = false;
	static bool first = true;
	
	char *lstype = config.LSTYPE;

	bool wait = false;
	for (int i=0; ;i++) {
		char c = lstype[i];
		if (c==0) break;
		if (c=='w') {
			wait = true;
			continue;
		}
		switch (c) {
			case 'v': RUN_V = true; break;
			case 'u': RUN_U = true; break;
			case 'q': RUN_Q = true; break;
			default: fprintf (stderr, "WARNING: invalid local search parameter (%c).\n", c);
		}
		wait = false;
	}
    if (maxrounds < 0) maxrounds = 999; //999; //large number

	int n = g.VertexCount();
	EdgeCost oldcost = solution.GetCost();
	int rounds = 0;
	int i=0;
	for (i=0; i<maxrounds; i++) {
		if(GetTime()-begin_time>config.TIME_LIMIT)
			break;
		rounds ++;
		if (RUN_V) {
			LSVertexInsertion::VertexInsertion(g, solution, n);
			MSTPrune(g, solution);
			if (verbose) fprintf (stderr, " v%d ", solution.GetCost());
		}
		//fprintf (stderr, "Starting key vertex elimination!\n");
		//fflush (stderr);
		
		if (RUN_Q) {
			LSKeyPath::KeyVertexElimination(g, solution);
			//fprintf (stderr, "Ending key vertex elimination!\n");
			if (verbose) fprintf (stderr, " q%d ", solution.GetCost());
		} 
		
		if (RUN_U) {
			LSVertexElimination::VertexElimination(g, solution);
			if (verbose) fprintf (stderr, " u%d ", solution.GetCost());
		}

		EdgeCost newcost = solution.GetCost();
		if (newcost > oldcost - EDGE_COST_PRECISION) break; //stop if result did not improve
		executionLogPtr.AddSolution(solution);
		oldcost = newcost;
	}
	first = false;
}

void RunLocalSearch (Graph &g, SteinerSolution &solution, int maxrounds, SteinerConfig &config) {
	bool RUN_Q = false;
	bool RUN_V = false;
	bool RUN_U = false;
	bool RUN_K = false;
	bool RESTRICT_K = false; //
	bool verbose = false;
	static bool first = true;
	
	char *lstype = config.LSTYPE;

	bool wait = false;
	for (int i=0; ;i++) {
		char c = lstype[i];
		if (c==0) break;
		if (c=='w') {
			wait = true;
			continue;
		}
		switch (c) {
			case 'v': RUN_V = true; break;
			case 'u': RUN_U = true; break;
			case 'q': RUN_Q = true; break;
			default: fprintf (stderr, "WARNING: invalid local search parameter (%c).\n", c);
		}
		wait = false;
	}
    if (maxrounds < 0) maxrounds = 999; //999; //large number

	int n = g.VertexCount();
	EdgeCost oldcost = solution.GetCost();
	int rounds = 0;
	int i=0;
	for (i=0; i<maxrounds; i++) {
		if(GetTime()-begin_time>config.TIME_LIMIT)
			break;
		rounds ++;
		if (RUN_V) {
			LSVertexInsertion::VertexInsertion(g, solution, n);
			MSTPrune(g, solution);
			if (verbose) fprintf (stderr, " v%d ", solution.GetCost());
		}
		//fprintf (stderr, "Starting key vertex elimination!\n");
		//fflush (stderr);
		
		if (RUN_Q) {
			LSKeyPath::KeyVertexElimination(g, solution);
			//fprintf (stderr, "Ending key vertex elimination!\n");
			if (verbose) fprintf (stderr, " q%d ", solution.GetCost());
		} 
		
		if (RUN_U) {
			LSVertexElimination::VertexElimination(g, solution);
			if (verbose) fprintf (stderr, " u%d ", solution.GetCost());
		}

		EdgeCost newcost = solution.GetCost();
		if (newcost > oldcost - EDGE_COST_PRECISION) break; //stop if result did not improve
		oldcost = newcost;
	}
	first = false;
}




static void DecayLocalSearch (SteinerSolution &solution, vector<EdgeCost> &current,vector<EdgeCost> &original, int decaysteps, double exponent, SteinerConfig &config) {
	Graph &g = *solution.g;
	int m = g.EdgeCount();
	const bool verbose = false;

	int decaytogo = decaysteps;
	
	for (;;) {
		EdgeCost prevcost = solution.GetCost();
		//cout<<endl<<decaytogo<<" *****  "<<prevcost<<endl;
		// run one round of local seach
		MSTPrune(g,solution);

		RunLocalSearch(g, solution, 1, config);
		

		EdgeCost newcost = solution.GetCost();
		
		if (decaytogo == 0) break;
		decaytogo --;

		for (int e=1; e<=m; e++) {
			//current[e] = original[e] + (current[e] - original[e]) * exponent; //NOT REALLY AN EXPONENT
			current[e] = original[e] + (current[e] - original[e]) * exponent; //NOT REALLY AN EXPONENT
		}
		g.ApplyCosts(current);
		solution.UpdateCost();
	}

	g.ApplyCosts(original); //restore original edge costs
	solution.UpdateCost();  //recost the existing solution
	EdgeCost before = solution.GetCost(); // this is the original cost
	
	// run local search on current solution using original costs 
	MSTPrune(g,solution);
	RunLocalSearch(g, solution, -1, config);		
	EdgeCost after = solution.GetCost();
	if (verbose) {
		fprintf (stderr, " %.0f", after);
	}
}

static void GenerateRandomizedSolution(SteinerSolution &solution, int root, vector<EdgeCost> &pertcost, SteinerConfig &config) 
{
	Graph &g = *solution.g;
	int m = g.EdgeCount();
			
	vector<EdgeCost> original (m+1);
	//remember original costs
	g.RetrieveCosts(original);

	g.ApplyCosts(pertcost);
	ConstructiveAlgorithms::SPH (g, solution, NULL, root);

	int LS_PERT_ROUNDS= max(999,g.VertexCount()); 
	double LS_PERT_EXPONENT = 1.0; //no decay
	
	LS_PERT_ROUNDS = config.LS_PERT_ROUNDS;
	LS_PERT_EXPONENT = config.LS_PERT_EXPONENT; 
	
	DecayLocalSearch(solution, pertcost, original,  LS_PERT_ROUNDS, LS_PERT_EXPONENT, config);
}





static void CombineSolutions(SteinerSolution &target, SteinerSolution &sa, SteinerSolution &sb, SteinerConfig &config, ExecutionLog &executionLogPtr) {
	Graph &g = *sa.g;
	int n = g.VertexCount();
	int m = g.EdgeCount();

	const bool verbose = false;

	vector<EdgeCost> pertcost(m+1,-1);

	for (int e = 1; e <= m; e++) {
		int tcount = 0;
		if (sa.Contains(e)) tcount++;
		if (sb.Contains(e)) tcount++;

		int mult = 1;
		if (tcount == 0) mult = 1000; //random.GetInteger(200,300); //edge in neither solution: very expensive
		else if (tcount == 1) mult = rand()%400+100;  //split edge: intermediate cost
		else { mult = 1; } //random.GetInteger(100,200); } //edge in both: keep it
		pertcost[e] = g.GetCost(e) * mult; // g.GetCost(a);
	}

	int root = Basics::PickRandomTerminal(g);
	ConstructiveAlgorithms::SPH (g, target, &pertcost[0], root);
	MSTPrune(g,target);
	RunLocalSearch(g, target, -1, config, executionLogPtr);

	if (verbose) fprintf (stderr, "%d x %d -> %d\n", sa.GetCost(), sb.GetCost(), target.GetCost());

}


static void CascadedCombination(SteinerSolution &solution, SteinerSolution &combsol, SolutionPool &elite, int maxfail,  SteinerConfig &config, ExecutionLog &executionLogPtr) {
	if (maxfail < 0) maxfail = config.MAX_COMB_FAIL;
	int failures_to_go = maxfail;
	const bool verbose = false;

	if (verbose) fprintf (stderr, "%d->", solution.GetCost());
	while (failures_to_go > 0) {

		if(GetTime()-begin_time>config.TIME_LIMIT)
			break;

		int rf=rand()%(elite.Count())+1;
		SteinerSolution *refsol = elite.GetReference(rf);
		CombineSolutions(combsol, solution, *refsol, config, executionLogPtr);
		if (!combsol.IsBetter(&solution)) {
			failures_to_go --;
		} else {
			solution.CopyFrom(&combsol);
		}
	}
	if (verbose) fprintf (stderr, "%d\n", solution.GetCost());
	elite.Add(&solution);
}


static void FlexibleMultistart(SteinerSolution &solution, SolutionPool &elite, int COMBINATION_THRESHOLD, SteinerConfig &config, ExecutionLog &executionLogPtr,vector<int> &edgeappear) {
	Graph &g = *solution.g;
	bool USE_PERTURBATION=true;
	
	int n = g.VertexCount();
	int m = g.EdgeCount();

	SteinerSolution bestsol(&g); //best solution found so far
	SteinerSolution combsol(&g); //combined solution
	bestsol.CopyFrom(&solution);
	

	EdgeCost bestcost = INFINITE_COST;
	if (elite.GetCount() > 0) {
		int p = elite.FindBestPosition();
		bestsol.CopyFrom(elite.GetReference(p));
		bestcost = bestsol.GetCost();
	}

	vector<EdgeCost> pertcost (m+1,-1);
	
	int root = rand()%n+1; //PickRandomTerminal(g, random);
	
	if (USE_PERTURBATION) {
		if (config.ga==0) {
			PerturbationTools::PartitionPerturbation(g,pertcost,edgeappear,config);
		} else {
			PerturbationTools::InitPerturbation(g, pertcost, config);
		}
	}

	GenerateRandomizedSolution (solution, root, pertcost, config);

	executionLogPtr.AddSolution(solution);
	
	elite.Add(&solution);

	if (elite.Count() == COMBINATION_THRESHOLD) {
		CascadedCombination(solution, combsol, elite, -1, config, executionLogPtr);
	}

	EdgeCost solcost = solution.GetCost();
	if (solcost < bestcost) {
		executionLogPtr.AddSolution(solution);
		bestcost = solcost;
		bestsol.CopyFrom(&solution);
	}

	solution.CopyFrom(&bestsol);
}

EdgeCost RunMultistart(Graph &g, SteinerConfig &config, ExecutionLog &executionLogPtr, SteinerSolution *outSolution, vector<int> &edgeappear) {
	double mstime = 0;
    double eps=0.000001;
	int msit=config.msit;

	EdgeCost bestfound = INFINITE_COST; 
    bool USE_PERTURBATION=false;
    int m=g.EdgeCount();
	long long doneit = 0;
	SolutionPool elite (100);
	int COMBINATION_THRESHOLD=1;
	SteinerSolution solution(&g);

	if (g.EdgeCount()>0 && g.TerminalCount()>1){
		while(1)
		{
			doneit++;
			
			COMBINATION_THRESHOLD=min((int)sqrt(doneit)+1,100);
			elite.curcap=COMBINATION_THRESHOLD;
			//cout<<"Time: "<<GetTime()-begin_time <<"    doneit: "<<doneit<<"    thresholds: "<<COMBINATION_THRESHOLD<<endl;
            FlexibleMultistart (solution, elite, COMBINATION_THRESHOLD, config, executionLogPtr, edgeappear);

			if (solution.GetCost()+eps < bestfound) {
				fprintf(stderr, "Time: %.2f  FOUND BETTER SOLUTION. IMPROVING COST FROM %.12f TO %.12f.\n",GetTime()-begin_time,bestfound, solution.GetCost());
				bestfound = solution.GetCost();
				if (outSolution != nullptr && (outSolution->GetCost() > solution.GetCost() || outSolution->EdgeCount() == 0)) {
					outSolution->CopyFrom(&solution);
				}
			}
			if(GetTime()-begin_time > config.TIME_LIMIT)
				break;
		}
	}
	else bestfound=g.GetFixedCost();

	fprintf (stderr, "Solution is %.12f.\n", (double)bestfound);
	return bestfound;
}




/*-------------------------------------for Partition Search----------------------------------------------*/
static void FlexibleMultistartforPartition(SteinerSolution &solution, SolutionPool &elite,int COMBINATION_THRESHOLD, SteinerConfig &config, ExecutionLog &executionLogPtr) {
	Graph &g = *solution.g;
	bool USE_PERTURBATION=true;
	
	int n = g.VertexCount();
	int m = g.EdgeCount();

	SteinerSolution bestsol(&g); //best solution found so far
	SteinerSolution combsol(&g); //combined solution
	bestsol.CopyFrom(&solution);
	

	EdgeCost bestcost = INFINITE_COST;
	if (elite.GetCount() > 0) {
		int p = elite.FindBestPosition();
		bestsol.CopyFrom(elite.GetReference(p));
		bestcost = bestsol.GetCost();
	}

	vector<EdgeCost> pertcost (m+1,-1);
	
	int root = rand()%n+1; //PickRandomTerminal(g, random);
	
	if (USE_PERTURBATION) {
		PerturbationTools::InitPerturbation(g, pertcost, config);
	}

	GenerateRandomizedSolution (solution, root, pertcost, config);

	executionLogPtr.AddSolution(solution);
	
	elite.Add(&solution);

	if (elite.Count() == COMBINATION_THRESHOLD) {
		CascadedCombination(solution, combsol, elite, -1, config, executionLogPtr);
	}

	EdgeCost solcost = solution.GetCost();
	if (solcost < bestcost) {
		executionLogPtr.AddSolution(solution);
		bestcost = solcost;
		bestsol.CopyFrom(&solution);
	}

	solution.CopyFrom(&bestsol);
}

EdgeCost RunMultistartforPartition(Graph &g, SteinerConfig &config, ExecutionLog &executionLogPtr, SteinerSolution *outSolution) {
	double mstime = 0;
    double eps=0.000001;
	int msit=config.msit;

	EdgeCost bestfound = INFINITE_COST; 
    bool USE_PERTURBATION=false;
    int m=g.EdgeCount();
	long long doneit = 0;
	SolutionPool elite (100);
	int COMBINATION_THRESHOLD=1;
	SteinerSolution solution(&g);

	if (g.EdgeCount()>0 && g.TerminalCount()>1){
		while(1)
		{
			doneit++;
			
			COMBINATION_THRESHOLD=min((int)sqrt(doneit)+1,100);
			elite.curcap=COMBINATION_THRESHOLD;
			//cout<<"Time: "<<GetTime()-begin_time <<"    doneit: "<<doneit<<"    thresholds: "<<COMBINATION_THRESHOLD<<endl;
            FlexibleMultistartforPartition (solution, elite, COMBINATION_THRESHOLD, config, executionLogPtr);

			if (solution.GetCost()+eps < bestfound) {
				fprintf(stderr, "Time: %.2f  FOUND BETTER SOLUTION. IMPROVING COST FROM %.12f TO %.12f.\n",GetTime()-begin_time,bestfound, solution.GetCost());
				bestfound = solution.GetCost();
				if (outSolution != nullptr && (outSolution->GetCost() > solution.GetCost() || outSolution->EdgeCount() == 0)) {
					outSolution->CopyFrom(&solution);
				}
			}
			if(GetTime()-begin_time > config.TIME_LIMIT)
				break;
		}
	}
	else bestfound=g.GetFixedCost();

	fprintf (stderr, "Solution is %.12f.\n", (double)bestfound);
	return bestfound;
}