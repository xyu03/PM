#include "timer.h"
#include "config.h"
#include "graph.h"
#include "solution.h"
#include "execution_log.h"
#include "stack.h"
#include "voronoi.h"
#include "uset.h"
#include "uf.h"
#include "basic_function.h"
#include "preprocess.h"
#include "constructive.h"
#include "LSVertexInsertion.h"
#include "LSKeyPath.h"
#include "LSVertexElimination.h"
#include "solver.h"
#include "connect_mincut.h"
#include "connect_shortestpath.h"
#include "connect_random.h"
#include "partition.h"
using namespace std;

void SolveInstance(int argc,char **argv)
{
    Graph g;
	g.ReadSTP(argv[1]);
    char *filename = argv[1];
    SteinerConfig config;

    config.ReadParameter(argc,argv);
    GetBeginTime();
    
    srand(config.seed+17);
    if (config.PREPROCESS) {
		fprintf (stderr, "Should be preprocessing.\n");
		double begin_prep_time=GetTime();
		
		Preprocessing::RunPreprocessing(g, !config.SAFE_PREPROCESS);
		
		fprintf (stdout, "preptime %.6f\n",GetTime()-begin_prep_time);
	}
	
	SteinerSolution bestSolution(&g);
	PartitionSearch(g, config, filename);
}   

int main(int argc,char **argv)
{
    SolveInstance(argc,argv);
    system("pause");
    return 0;
}
