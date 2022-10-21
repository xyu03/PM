vector<int> fathersp;
vector<int> sizsp;
vector<int> terminalnodessp;

struct nodebridgesp{
    int eid;
    double dis;
};
bool cmpbridgessp(const nodebridgesp &a,nodebridgesp &b)
{
    return a.dis<b.dis;
}

vector<nodebridgesp>bridgeedgessp;
int GetFatherforPartitionSP(int v)
{
    return v==fathersp[v]?v:fathersp[v]=GetFatherforPartitionSP(fathersp[v]);
}

void PartitionPreProcessSP(Graph &g,SteinerConfig &config,VoronoiData &vordata, UniverseSet &baselist, BinaryHeap<EdgeCost> &binheap)
{
    int n=g.VertexCount();
    int m=g.EdgeCount();
    int t=g.TerminalCount();
    config.partnumbase=min(config.partnumbase,t);

    if(config.partnum==0)
        config.partnum=t/config.partnumbase;

    Basics::ComputeVoronoi(g,vordata,baselist,binheap,NULL); 

    fathersp.resize(n+1,0); 
    sizsp.resize(n+1,0);

    terminalnodessp.clear();
    for(int v=1;v<=n;v++){
        if(g.IsTerminal(v))
            terminalnodessp.push_back(v);
    }

    bridgeedgessp.clear();
    nodebridgesp nb;
    for(int i=1;i<=m;i++)
    {
        int u,v;
        g.GetEndpoints(i,u,v);
        if(vordata.GetBase(u)!=vordata.GetBase(v)){
            nb.eid=i;
            nb.dis=vordata.GetDistance(u)+vordata.GetDistance(v)+g.GetCost(i);
            bridgeedgessp.push_back(nb);
        }
    }
}



int PartitionSP(Graph &g,SteinerConfig &config,VoronoiData &vordata)
{
    int n=g.VertexCount();
    int m=g.EdgeCount();
    int t=g.TerminalCount();

    int k=config.partnum;
    int d=1;/////
    int curk=t; 
    int bridgenum=bridgeedgessp.size();

    
    for(int v=0;v<=n;v++){
        fathersp[v]=0;
        sizsp[v]=0;
    }

    for(int v=1;v<=n;v++){
        int fv=vordata.GetBase(v);
        fathersp[v]=fv;
        sizsp[fv]++;
    }
    for(int i=t-1;i>=0;i--){
        swap(terminalnodessp[i],terminalnodessp[rand()%(i+1)]);
    }


    sort(bridgeedgessp.begin(),bridgeedgessp.end(),cmpbridgessp);
    int fbase=2;
    for(int f=t;f>=1;f=f/fbase)
    {
        for(int i=0;i<bridgenum;i++)
        {
            int u,v;
            int e=bridgeedgessp[i].eid;
            g.GetEndpoints(e,u,v);
            int fu=GetFatherforPartitionSP(u);
            int fv=GetFatherforPartitionSP(v);
            if(fu==fv) continue;
            if(sizsp[fu]+sizsp[fv]<=n/f*d)
            {
                fathersp[fu]=fv;
                sizsp[fv]+=sizsp[fu];
                curk--;
                if(curk==k) 
                    break;
            }
        }
        if(curk==k) break;
    }

    for(int i=1;i<=g.VertexCount();i++)
    {
        fathersp[i]=GetFatherforPartitionSP(i);
        g.SetBelongPartition(i,fathersp[i]);
    }


    bool verbox=true;
    if(verbox)
    {
        printf("-----------------initial partition info---------------------------\n");
        int num=0;
        int it=0;
        for(int i=0;i<t;i++)
        {
            int v=terminalnodessp[i];
            int fv=GetFatherforPartitionSP(v);
            if(fv!=v) continue;
            int tcount=0;
            for(int j=0;j<t;j++)
            {
                int vv=terminalnodessp[j];
                vv=GetFatherforPartitionSP(vv);
                if(vv==v)
                    tcount++;
            }
            printf("%d partition %d   %d   %d\n",++it,v,sizsp[v],tcount);
            num+=sizsp[v];
        }
        printf("partition num:  %d\n",curk);
        printf("total vertex num:  %d\n",num);
        printf("-----------------initial partition info---------------------------\n");
    }
    return curk;
}

void PMSearchSP(Graph &g,SteinerConfig &config,int partnum,vector<int>&zeroedge)
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
    for(int i=0;i<terminalnodessp.size();i++)
    {
        int v=terminalnodessp[i];
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


    vector<nodebridgesp>bridgeedgesspML;
    bridgeedgesspML.clear();
    for(int i=0;i<bridgeedgessp.size();i++)
    {
        int e=bridgeedgessp[i].eid;
        int u,v,fu,fv;
        g.GetEndpoints(e,u,v);
        fu=GetFatherforPartitionSP(u);
        fv=GetFatherforPartitionSP(v);
        if(fu!=fv){
            bridgeedgesspML.push_back(bridgeedgessp[i]);
        }
    }
    int bridgeedgesspMLnum=bridgeedgesspML.size();
    sort(bridgeedgesspML.begin(),bridgeedgesspML.end(),cmpbridgessp);

    cout<<"bridgeML  num:  "<<bridgeedgesspMLnum<<endl;

    vector<int>zeroedgesub;
    zeroedgesub.resize(m,0);

    int curk=partnum;
    vector<bool>ifcombined;
    ifcombined.resize(n,false);
    int fbase=config.fbase;
    for(int f=k-1;f>=1;f=f/fbase)
    {
        for(int i=0;i<bridgeedgesspMLnum;i++)
        {
            int u,v;
            int e=bridgeedgesspML[i].eid;
            g.GetEndpoints(e,u,v);
            int fu=GetFatherforPartitionSP(u);
            int fv=GetFatherforPartitionSP(v);
            if(fu==fv) continue;
            if(sizsp[fu]+sizsp[fv]<=n/f*d)
            {
                fathersp[fu]=fv;
                ifcombined[fv]=true;
                sizsp[fv]+=sizsp[fu];
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
            fathersp[i]=GetFatherforPartitionSP(i);
            g.SetBelongPartition(i,fathersp[i]);
        }

        int partitionindex=0;
        config.TIME_LIMIT=GetTime()-begin_time;
        for(int i=0;i<terminalnodessp.size();i++)
        {
            int v=terminalnodessp[i];
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
    }
}