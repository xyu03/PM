
vector<int> fatherrd;
vector<int> sizrd;
vector<int> terminalnodesrd;

struct nodebridgerd{
    int eid;
    long long dis;
};
bool cmpbridges(const nodebridgerd &a,nodebridgerd &b)
{
    return a.dis<b.dis;
}

vector<nodebridgerd>bridgeedgesrd;
vector<nodebridgerd>bridgeedgesrd1;
int GetFatherforPartitionRandom(int v)
{
    return v==fatherrd[v]?v:fatherrd[v]=GetFatherforPartitionRandom(fatherrd[v]);
}

void PartitionPreProcessRandom(Graph &g,SteinerConfig &config,VoronoiData &vordata, UniverseSet &baselist, BinaryHeap<EdgeCost> &binheap)
{
    int n=g.VertexCount();
    int m=g.EdgeCount();
    int t=g.TerminalCount();
    config.partnumbase=min(config.partnumbase,t);

    if(config.partnum==0)
        config.partnum=t/config.partnumbase;

    Basics::ComputeVoronoi(g,vordata,baselist,binheap,NULL); 

    fatherrd.resize(n+1,0); 
    sizrd.resize(n+1,0);

    terminalnodesrd.clear();
    for(int v=1;v<=n;v++){
        if(g.IsTerminal(v))
            terminalnodesrd.push_back(v);
    }

    bridgeedgesrd1.clear();
    nodebridgerd nb;
    for(int i=1;i<=m;i++)
    {
        int u,v;
        g.GetEndpoints(i,u,v);
        if(vordata.GetBase(u)!=vordata.GetBase(v)){
            nb.eid=i;
            int fu=vordata.GetBase(u);
            int fv=vordata.GetBase(v);
            if(fu>fv) swap(fu,fv);
            nb.dis=10000000ll*fu+fv;
            bridgeedgesrd1.push_back(nb);
        }
    }
    sort(bridgeedgesrd1.begin(),bridgeedgesrd1.end(),cmpbridges);
    bridgeedgesrd.clear();
    if(bridgeedgesrd1.size())
        bridgeedgesrd.push_back(bridgeedgesrd1[0]);
    for(int i=1;i<bridgeedgesrd1.size();i++)
    {
        if(bridgeedgesrd1[i].dis!=bridgeedgesrd1[i-1].dis)
            bridgeedgesrd.push_back(bridgeedgesrd1[i]);
    }

}



int PartitionRandom(Graph &g,SteinerConfig &config,VoronoiData &vordata)
{
    int n=g.VertexCount();
    int m=g.EdgeCount();
    int t=g.TerminalCount();

    int k=config.partnum;
    int d=1;/////
    int curk=t; 
    int bridgenum=bridgeedgesrd.size();

    
    for(int v=0;v<=n;v++){
        fatherrd[v]=0;
        sizrd[v]=0;
    }

    for(int v=1;v<=n;v++){
        int fv=vordata.GetBase(v);
        fatherrd[v]=fv;
        sizrd[fv]++;
    }
    for(int i=t-1;i>=0;i--){
        swap(terminalnodesrd[i],terminalnodesrd[rand()%(i+1)]);
    }

    for(int i=bridgenum-1;i>=0;i--){
        swap(bridgeedgesrd[i],bridgeedgesrd[rand()%(i+1)]);
    }
    
    int fbase=config.fbase;
    for(int f=t;f>=1;f=f/fbase)
    {
        for(int i=0;i<bridgenum;i++)
        {
            int u,v;
            int e=bridgeedgesrd[i].eid;
            g.GetEndpoints(e,u,v);
            int fu=GetFatherforPartitionRandom(u);
            int fv=GetFatherforPartitionRandom(v);
            if(fu==fv) continue;
            if(sizrd[fu]+sizrd[fv]<=n/f*d)
            {
                fatherrd[fu]=fv;
                sizrd[fv]+=sizrd[fu];
                curk--;
                if(curk==k) 
                    break;
            }
        }
        if(curk==k) break;

        nodebridgerd nb;
        bridgeedgesrd1.clear();
        for(int i=0;i<bridgenum;i++){
            int u,v;
            int e=bridgeedgesrd[i].eid;
            g.GetEndpoints(e,u,v);
            int fu=GetFatherforPartitionRandom(u);
            int fv=GetFatherforPartitionRandom(v);
            if(fu==fv) continue;
            
            if(fu>fv) swap(fu,fv);
            nb.dis=10000000ll*fu+fv;
            nb.eid=e;
            bridgeedgesrd1.push_back(nb);
        }
        sort(bridgeedgesrd1.begin(),bridgeedgesrd1.end(),cmpbridges);

        bridgeedgesrd.clear();
        if(bridgeedgesrd1.size())
            bridgeedgesrd.push_back(bridgeedgesrd1[0]);
        for(int i=1;i<bridgeedgesrd1.size();i++)
        {
            if(bridgeedgesrd1[i].dis!=bridgeedgesrd1[i-1].dis)
                bridgeedgesrd.push_back(bridgeedgesrd1[i]);
        }

        bridgenum=bridgeedgesrd.size();

        for(int i=bridgenum-1;i>=0;i--){
            swap(bridgeedgesrd[i],bridgeedgesrd[rand()%(i+1)]);
        }
    }

    for(int i=1;i<=g.VertexCount();i++)
    {
        fatherrd[i]=GetFatherforPartitionRandom(i);
        g.SetBelongPartition(i,fatherrd[i]);
    }


    bool verbox=true;
    if(verbox)
    {
        printf("-----------------initial partition info---------------------------\n");
        int num=0;
        int it=0;
        for(int i=0;i<t;i++)
        {
            int v=terminalnodesrd[i];
            int fv=GetFatherforPartitionRandom(v);
            if(fv!=v) continue;
            int tcount=0;
            for(int j=0;j<t;j++)
            {
                int vv=terminalnodesrd[j];
                vv=GetFatherforPartitionRandom(vv);
                if(vv==v)
                    tcount++;
            }
            printf("%d partition %d   %d   %d\n",++it,v,sizrd[v],tcount);
            num+=sizrd[v];
        }
        printf("partition num:  %d\n",curk);
        printf("total vertex num:  %d\n",num);
        printf("-----------------initial partition info---------------------------\n");
    }
    return curk;
}

void PMSearchRandom(Graph &g,SteinerConfig &config,int partnum,vector<int>&zeroedge)
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
    for(int i=0;i<terminalnodesrd.size();i++)
    {
        int v=terminalnodesrd[i];
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


    vector<nodebridgerd>bridgeedgesrdML;
    bridgeedgesrdML.clear();
    for(int i=0;i<bridgeedgesrd.size();i++)
    {
        int e=bridgeedgesrd[i].eid;
        int u,v,fu,fv;
        g.GetEndpoints(e,u,v);
        fu=GetFatherforPartitionRandom(u);
        fv=GetFatherforPartitionRandom(v);
        if(fu!=fv){
            bridgeedgesrdML.push_back(bridgeedgesrd[i]);
        }
    }
    int bridgeedgesrdMLnum=bridgeedgesrdML.size();
    for(int i=bridgeedgesrdMLnum-1;i>=0;i--){
        swap(bridgeedgesrdML[i],bridgeedgesrdML[rand()%(i+1)]);
    }

    cout<<"bridgeML  num:  "<<bridgeedgesrdMLnum<<endl;

    vector<int>zeroedgesub;
    zeroedgesub.resize(m,0);

    int curk=partnum;
    vector<bool>ifcombined;
    ifcombined.resize(n,false);
    int fbase=config.fbase;
    for(int f=k-1;f>=1;f=f/fbase)
    {
        for(int i=0;i<bridgeedgesrdMLnum;i++)
        {
            int u,v;
            int e=bridgeedgesrdML[i].eid;
            g.GetEndpoints(e,u,v);
            int fu=GetFatherforPartitionRandom(u);
            int fv=GetFatherforPartitionRandom(v);
            if(fu==fv) continue;
            if(sizrd[fu]+sizrd[fv]<=n/f*d)
            {
                fatherrd[fu]=fv;
                ifcombined[fv]=true;
                sizrd[fv]+=sizrd[fu];
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
            fatherrd[i]=GetFatherforPartitionRandom(i);
            g.SetBelongPartition(i,fatherrd[i]);
        }

        int partitionindex=0;
        config.TIME_LIMIT=GetTime()-begin_time;
        for(int i=0;i<terminalnodesrd.size();i++)
        {
            int v=terminalnodesrd[i];
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
        nodebridgerd nb;

        bridgeedgesrd1.clear();
        for(int i=0;i<bridgeedgesrdMLnum;i++){
            int u,v;
            int e=bridgeedgesrdML[i].eid;
            g.GetEndpoints(e,u,v);
            int fu=GetFatherforPartitionRandom(u);
            int fv=GetFatherforPartitionRandom(v);
            if(fu==fv) continue;
            
            if(fu>fv) swap(fu,fv);
            nb.dis=10000000ll*fu+fv;
            nb.eid=e;
            bridgeedgesrd1.push_back(nb);
        }
        sort(bridgeedgesrd1.begin(),bridgeedgesrd1.end(),cmpbridges);
        
        bridgeedgesrdML.clear();
        if(bridgeedgesrd1.size())
            bridgeedgesrdML.push_back(bridgeedgesrd1[0]);
        for(int i=1;i<bridgeedgesrd1.size();i++)
        {
            if(bridgeedgesrd1[i].dis!=bridgeedgesrd1[i-1].dis)
                bridgeedgesrdML.push_back(bridgeedgesrd1[i]);
        }

        bridgeedgesrdMLnum=bridgeedgesrdML.size();

        for(int i=bridgeedgesrdMLnum-1;i>=0;i--){
            swap(bridgeedgesrdML[i],bridgeedgesrdML[rand()%(i+1)]);
        }
    }
}
