vector<int> fathermincut;
vector<int> sizmincut;
vector<int> terminalnodesmincut;
vector<int>bridgeedgesmincut;
int GetFatherforPartitionMincut(int v)
{
    return v==fathermincut[v]?v:fathermincut[v]=GetFatherforPartitionMincut(fathermincut[v]);
}

void PartitionPreProcessMincut(Graph &g,SteinerConfig &config,VoronoiData &vordata, UniverseSet &baselist, BinaryHeap<EdgeCost> &binheap)
{
    int n=g.VertexCount();
    int m=g.EdgeCount();
    int t=g.TerminalCount();
    config.partnumbase=min(config.partnumbase,t);

    if(config.partnum==0)
        config.partnum=t/config.partnumbase;

    Basics::ComputeVoronoi(g,vordata,baselist,binheap,NULL); 

    fathermincut.resize(n+1,0); 
    sizmincut.resize(n+1,0);

    terminalnodesmincut.clear();
    for(int v=1;v<=n;v++){
        if(g.IsTerminal(v))
            terminalnodesmincut.push_back(v);
    }

    bridgeedgesmincut.clear();
    for(int i=1;i<=m;i++)
    {
        int u,v;
        g.GetEndpoints(i,u,v);
        if(vordata.GetBase(u)!=vordata.GetBase(v)){
            bridgeedgesmincut.push_back(i);
        }
    }
}

int PartitionMincut(Graph &g,SteinerConfig &config,VoronoiData &vordata)
{
    int n=g.VertexCount();
    int m=g.EdgeCount();
    int t=g.TerminalCount();

    int k=config.partnum;
    int d=1;/////
    int curk=t; 
    int bridgenum=bridgeedgesmincut.size();

    
    for(int v=0;v<=n;v++){
        fathermincut[v]=0;
        sizmincut[v]=0;
    }

    for(int v=1;v<=n;v++){
        int fv=vordata.GetBase(v);
        fathermincut[v]=fv;
        sizmincut[fv]++;
    }
    for(int i=t-1;i>=0;i--){
        swap(terminalnodesmincut[i],terminalnodesmincut[rand()%(i+1)]);
    }


    for(int i=bridgenum-1;i>=0;i--)
    {
        swap(bridgeedgesmincut[i],bridgeedgesmincut[rand()%(i+1)]);
    }
    int fbase=config.fbase;
    for(int f=t;f>=1;f=f/fbase)
    {
        for(int i=0;i<bridgenum;i++)
        {
            int u,v;
            int e=bridgeedgesmincut[i];
            g.GetEndpoints(e,u,v);
            int fu=GetFatherforPartitionMincut(u);
            int fv=GetFatherforPartitionMincut(v);
            if(fu==fv) continue;
            if(sizmincut[fu]+sizmincut[fv]<=n/f*d)
            {
                fathermincut[fu]=fv;
                sizmincut[fv]+=sizmincut[fu];
                curk--;
                if(curk==k) 
                    break;
            }
        }
        if(curk==k) break;

        for(int i=bridgenum-1;i>=0;i--){
            swap(bridgeedgesmincut[i],bridgeedgesmincut[rand()%(i+1)]);
        }
    }

    for(int i=1;i<=g.VertexCount();i++)
    {
        fathermincut[i]=GetFatherforPartitionMincut(i);
        g.SetBelongPartition(i,fathermincut[i]);
    }


    bool verbox=true;
    if(verbox)
    {
        printf("-----------------initial partition info---------------------------\n");
        int num=0;
        int it=0;
        for(int i=0;i<t;i++)
        {
            int v=terminalnodesmincut[i];
            int fv=GetFatherforPartitionMincut(v);
            if(fv!=v) continue;
            int tcount=0;
            for(int j=0;j<t;j++)
            {
                int vv=terminalnodesmincut[j];
                vv=GetFatherforPartitionMincut(vv);
                if(vv==v)
                    tcount++;
            }
            printf("%d partition %d   %d   %d\n",++it,v,sizmincut[v],tcount);
            num+=sizmincut[v];
        }
        printf("partition num:  %d\n",curk);
        printf("total vertex num:  %d\n",num);
        printf("-----------------initial partition info---------------------------\n");
    }
    return curk;
}

void PMSearchMincut(Graph &g,SteinerConfig &config,int partnum,vector<int>&zeroedge)
{
    int n=g.VertexCount();
    int m=g.EdgeCount();
    int t=g.TerminalCount();
    int k=partnum;
    int d=1;

    double delttimelimit=config.delttime;
    //delttimelimit=1;
    EdgeCost sumofpartcost;
    int partitionindex=0;
    config.TIME_LIMIT=GetTime()-begin_time;
    for(int i=0;i<terminalnodesmincut.size();i++)
    {
        int v=terminalnodesmincut[i];
        int fv=g.GetBelong_to_Patition(v);
        if(fv!=v) continue;
        printf("---------------------begin to solve partition %d------------------------\n",++partitionindex);
        double beginsolvepartitiontime=GetTime();
        Graph g1;
        g1.ReadSubSTP(g,v);

        ExecutionLog executionLog(&g1, beginsolvepartitiontime, config.TIME_LIMIT);
        SteinerSolution bestSolution(&g1);

        config.TIME_LIMIT+=delttimelimit;
        EdgeCost partcost=RunMultistartforPartition(g1, config, executionLog, &bestSolution);
        sumofpartcost+=partcost;

        if (bestSolution.EdgeCount()!=0) {
            int v,w;
            for (size_t e = 1; e <= g1.EdgeCount(); ++e) 
            if (bestSolution.Contains(e))
            {
                zeroedge[oldedgeid[e]]=1;
            }
        }
        printf("-----------------------finished partition %d---------------------------\n\n",partitionindex);
    }


    vector<int>bridgeedgesmincutML;
    bridgeedgesmincutML.clear();
    for(int i=0;i<bridgeedgesmincut.size();i++)
    {
        int e=bridgeedgesmincut[i];
        int u,v,fu,fv;
        g.GetEndpoints(e,u,v);
        fu=GetFatherforPartitionMincut(u);
        fv=GetFatherforPartitionMincut(v);
        if(fu!=fv){
            bridgeedgesmincutML.push_back(e);
        }
    }
    int bridgeedgesmincutMLnum=bridgeedgesmincutML.size();
    for(int i=bridgeedgesmincutMLnum-1;i>=0;i--){
        swap(bridgeedgesmincutML[i],bridgeedgesmincutML[rand()%(i+1)]);
    }

    cout<<"bridgeML  num:  "<<bridgeedgesmincutMLnum<<endl;

    vector<int>zeroedgesub;
    zeroedgesub.resize(m,0);

    int curk=partnum;
    vector<bool>ifcombined;
    ifcombined.resize(n,false);
    int fbase=config.fbase;
    for(int f=k-1;f>=1;f=f/fbase)
    {
        for(int i=0;i<bridgeedgesmincutMLnum;i++)
        {
            int u,v;
            int e=bridgeedgesmincutML[i];
            g.GetEndpoints(e,u,v);
            int fu=GetFatherforPartitionMincut(u);
            int fv=GetFatherforPartitionMincut(v);
            if(fu==fv) continue;
            if(sizmincut[fu]+sizmincut[fv]<=n/f*d)
            {
                fathermincut[fu]=fv;
                ifcombined[fv]=true;
                sizmincut[fv]+=sizmincut[fu];
                curk--;
                if(curk==1) 
                    break;
            }
        }
        if(curk==1) break;

        ////SEARCH PART
        cout<<"hhhhhhhhhhhhhhhhhhhhhhhhhhhh"<<curk<<endl;

        for(int i=1;i<=g.VertexCount();i++)
        {
            fathermincut[i]=GetFatherforPartitionMincut(i);
            g.SetBelongPartition(i,fathermincut[i]);
        }

        int partitionindex=0;
        config.TIME_LIMIT=GetTime()-begin_time;
        for(int i=0;i<terminalnodesmincut.size();i++)
        {
            int v=terminalnodesmincut[i];
            int fv=g.GetBelong_to_Patition(v);
            if(fv!=v||!ifcombined[fv]) continue;
            ifcombined[fv]=false;
            
            printf("---------------------begin to solve partition %d------------------------\n",++partitionindex);
            double beginsolvepartitiontime=GetTime();
            Graph g1;
            g1.ReadSubSTP(g,v);
            for(int e=1;e<=g1.EdgeCount();e++)
            {
                zeroedgesub[e]=zeroedge[oldedgeid[e]];
            }

            config.TIME_LIMIT+=delttimelimit;
            ExecutionLog executionLog(&g1, beginsolvepartitiontime, config.TIME_LIMIT);
            SteinerSolution bestSolution(&g1);
            
            //EdgeCost partcost=RunMultistartforPartition(g1, config, executionLog, &bestSolution);
            EdgeCost partcost=RunMultistart(g1, config, executionLog, &bestSolution,zeroedgesub);
            sumofpartcost+=partcost;

            if (bestSolution.EdgeCount()!=0) {
                int v,w;
                for (size_t e = 1; e <= g1.EdgeCount(); ++e) 
                if (bestSolution.Contains(e)){
                    zeroedge[oldedgeid[e]]=1;
                }
                else zeroedge[oldedgeid[e]]=0;
            }
            printf("-----------------------finished partition %d---------------------------\n\n",partitionindex);
        }

        for(int i=bridgeedgesmincutMLnum-1;i>=0;i--){
           swap(bridgeedgesmincutML[i],bridgeedgesmincutML[rand()%(i+1)]);
        }
    }
}