

void PartitionSearch( Graph &g, SteinerConfig &config, char *filename)
{
    int partition_flag=(config.ga+1)%2;
    int n=g.VertexCount();
    int m=g.EdgeCount();
    int t=g.TerminalCount();
    double TIME_LIMIT=config.TIME_LIMIT;
    vector<int>zeroedge;
    zeroedge.resize(m+1,0);

    config.partnum=min(config.partnum,t/10);

    if(partition_flag==1)
    {
        printf("\nBEGIN PARTITION SEARCH !\n\n");
        UniverseSet baselist(n);
        for (int v=1; v<=n; v++) {
            if (g.IsTerminal(v)) baselist.Insert(v);
        }
        BinaryHeap <EdgeCost> binheap(n); // = new BinaryHeap<ArcCost> (n);
        VoronoiData vordata(n);// = new VoronoiData(n);

        if(config.connectmethod==1)///min-cut
        {
            PartitionPreProcessMincut(g,config,vordata,baselist,binheap);

            double begin_partition_time=GetTime();
            int partnum=PartitionMincut(g,config,vordata);

            PMSearchMincut(g,config,partnum,zeroedge);

            fprintf (stdout, "partition ptime %.6f\n",GetTime()-begin_partition_time);
        }

        if(config.connectmethod==2)///terminal shortest path
        {
            PartitionPreProcessSP(g,config,vordata,baselist,binheap);

            double begin_partition_time=GetTime();
            int partnum=PartitionSP(g,config,vordata);

            PMSearchSP(g,config,partnum,zeroedge);

            fprintf (stdout, "partition ptime %.6f\n",GetTime()-begin_partition_time);
        }

		if(config.connectmethod==3)///randomly
        {
            PartitionPreProcessRandom(g,config,vordata,baselist,binheap);

            double begin_partition_time=GetTime();
            int partnum=PartitionRandom(g,config,vordata);

            PMSearchRandom(g,config,partnum,zeroedge);

            fprintf (stdout, "partition ptime %.6f\n",GetTime()-begin_partition_time);
        }
        
        config.TIME_LIMIT=TIME_LIMIT;
        //printf("\nConnect partitions , BEGIN COMPLETE SEARCH !\n\n");
    }
    else printf("\nBEGIN LocalOptimize SEARCH !\n\n");

    ExecutionLog executionLog(&g, begin_time, config.TIME_LIMIT);
    SteinerSolution bestSolution(&g);
    
    EdgeCost answercost=RunMultistart(g, config, executionLog, &bestSolution, zeroedge);

    if (!config.LOG_FILENAME.empty()) {
		ofstream logFile(config.LOG_FILENAME.c_str());
		if (logFile.is_open()) {
			logFile << "SECTION Comment" << endl
				<< "Name \"" << filename << "\"" << endl
				<< "Problem \"SPG\"" << endl
				<< "Program \"PM\"" << endl
				<< "End" << endl
				<< endl
				<< "SECTION Solutions" << endl;
			for (size_t i = 0; i < executionLog.solCost.size(); ++i) {
				logFile << "Solution " << fixed << executionLog.solCost[i].second << " " << fixed <<setprecision(12)<< executionLog.solCost[i].first << endl;
			}
			logFile << "End" << endl
				<< endl
				<< "SECTION Run" << endl
				<< "Threads 1" << endl
				<< "Time " << GetTime()-begin_time << endl
				<< "Dual 0" << endl;
			logFile << "Primal " << fixed <<setprecision(12)<< answercost << endl;
			logFile << "End" << endl
				<< endl;

			if (!executionLog.solCost.empty()) {
				logFile << "SECTION Finalsolution" << endl;

				size_t numVertices = 0;
				int sum=0;
				for (int v = 1; v <= g.VertexCount(); ++v) {
					sum++;
					if (executionLog.bestSolution.GetDegree(v) > 0 || g.IsTerminal(v))
						++numVertices;
				}
				logFile << "Vertices " << numVertices << endl;
				for (int v = 1; v <= g.VertexCount(); ++v) {
					if (executionLog.bestSolution.GetDegree(v) > 0 || g.IsTerminal(v))
						logFile << "V " << v << endl;
				}
				logFile << "Edges " << executionLog.bestSolution.EdgeCount() << endl;
				for (size_t e = 1; e <= g.EdgeCount(); ++e) {
					if (executionLog.bestSolution.Contains(e))
						logFile << "E " << g.GetFirstEndpoint(e) << " " << g.GetSecondEndpoint(e) << endl;
				}
				logFile << "End" << endl;
			}
			logFile.close();
		}
	}
}