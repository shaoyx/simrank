#include "OIPSimRank.h"
void OIPSimRank::run(int qv, int k){
	if(isInit == false){
		isInit = true;
		initialize();
	}
	else{
	    char filepath[125];
	    sprintf(filepath, "dataset/%s/index/OIP.partial", graphName);
	    FILE* fp = fopen(filepath, "rb");
        if(fp == NULL){
            printf("Failed to open the %s file.\n", filepath);
        }
//        printf("reading allpair file: %s\n", filepath);
        for(int i = 0; i < maxVertexId; ++i){
	    	fread(srvalue[maxSteps&1][i], sizeof(double), maxVertexId, fp);
	    }
        fclose(fp);
	}

	vector<SimRankValue> res;
	res.resize(maxVertexId);
	for(int i = 0; i < maxVertexId; ++i){
		if(i > qv && srvalue[maxSteps & 1][qv][i] > 0.0){
			//res.push_back(SimRankValue(i, srvalue[maxSteps & 1][qv][i]));
			res[i].setVid(i);
            res[i].incValue(srvalue[maxSteps & 1][qv][i]);
		}
		else if(i < qv && srvalue[maxSteps & 1][i][qv] > 0.0){
			//res.push_back(SimRankValue(i, srvalue[maxSteps & 1][i][qv]));
			res[i].setVid(i);
            res[i].incValue(srvalue[maxSteps & 1][i][qv]);
		}
	}
   // printf("candidate size=%d\n", res.size());
	save(res, k);

}

void OIPSimRank::initialize(){
	//sort the outgoing neighbors for using intersecti
	#ifdef DEBUG
	printf("sort\n");
	#endif
	for(int i=0;i<maxVertexId;i++)
		sort(graphDst+graphSrc[i],graphDst+graphSrc[i+1]);
	#ifdef DEBUG
	printf("sort finished\n");
	#endif
    	generateWeightMatrix_MST();
	#ifdef DEBUG
	printf("mst finished\n");
	#endif

	int iter=0;
	double ** pSum = new double* [maxVertexId];
	for(int i=0;i<maxVertexId;i++){
		pSum[i]=new double [maxVertexId];//pSum[u][y] denotes pSum_Iu^y
	//	memset(pSum[i],0,sizeof(int)*maxVertexId);
	}
	#ifdef DEBUG
	printf("while loop start\n");
	#endif
	Time timer;
	timer.start();
	while(iter<maxSteps){
		for(int i=0;i<maxVertexId;i++)
			memset(pSum[i],0,sizeof(double)*maxVertexId);
		for(int i=0;i<maxVertexId;i++){
			//each edge in mst
			int v=mstSrc[i];
			int u=mstDst[i];
			int w=mstWeight[i];
			int scanCost= graphSrc[v+1]-graphSrc[v] + graphSrc[u+1]- graphSrc[u];
			int addCost = graphSrc[u+1]-graphSrc[u] - 1 - mstWeight[i];
			//double factor = 0.5;
			if(mstSrc[i]==maxVertexId || scanCost * decayFactor >addCost){// node maxVertexId represents node '#'
				int u=mstDst[i];//we compute parital[u][y] directly
				for(int y=0;y<maxVertexId;y++)
					for(int k=graphSrc[u];k<graphSrc[u+1];k++)
						pSum[u][y]+=srvalue[iter &1][y][graphDst[k]];
				OP(u,iter,pSum[u]);
			}
			else{//given mstSrc[i], compute mstDst[i]
			//	int v=mstSrc[i];
			//	int u=mstDst[i];
			//	int w=mstWeight[i];

				//use v to compute u
					for(int y=0;y<maxVertexId;y++){ 	
						pSum[u][y]=pSum[v][y];
					int vs=graphSrc[v],ve=graphSrc[v+1];
					int us=graphSrc[u],ue=graphSrc[u+1];
					while(true){
						if(us==ue){
							for(int vv=vs;vv<ve;vv++)
								pSum[u][y]-=srvalue[iter&1][y][graphDst[vv]];
							break;
						}
						if(vs==ve){
							for(int uu=us;uu<ue;uu++)
								pSum[u][y]+=srvalue[iter&1][y][graphDst[uu]];
							break;
						}
						if(graphDst[us]<graphDst[vs]){
							pSum[u][y]+=srvalue[iter&1][y][graphDst[us]];
							us++;
						}
						else if(graphDst[us]>graphDst[vs]){
							pSum[u][y]-=srvalue[iter&1][y][graphDst[vs]];
							vs++;
						}
						else{
							us++;vs++;
						}
					}	
				}	
				OP(u,iter,pSum[u]);
				
			}
	
		}

	iter++;
	}
	
	timer.stop();
	printf("time cost for while loop of OIP: %.5lf\n",timer.getElapsedTime());
	#ifdef DEBUG
	printf("while loop stopped\n");
	#endif
    char filepath[125];
    sprintf(filepath, "dataset/%s/index/OIP.partial", graphName);
    FILE* fp = fopen(filepath, "wb");
    for(int i = 0; i < maxVertexId; ++i){
    	fwrite(srvalue[maxSteps&1][i], sizeof(double), maxVertexId, fp);
    }
    fclose(fp);
	for( int i=0;i<maxVertexId;i++)
		delete [] pSum[i];
	delete [] pSum;
}
void OIPSimRank::generateWeightMatrix_MST(){
// output graphSrc,graphDst
	#ifdef DEBUG
	printf("output reversed Graph\n");
	for(int i=0;i<maxVertexId;i++)
		for(int j=graphSrc[i];j<graphSrc[i+1];j++)
		printf("mst edge(%d,%d)\n",i,graphDst[j]);
	#endif

	const static int maxWeight=100000;
	int **A=new int* [maxVertexId+1];
	for(int i=0;i<maxVertexId+1;i++){
		A[i]=new int[maxVertexId+1];//include node '#' as node N
		for(int j=0;j<maxVertexId+1;j++)
			A[i][j]=maxWeight;
	}
	//gennerate weight matrix store in A
	for(int i=0;i<maxVertexId;i++)
		for(int j=i+1;j<maxVertexId;j++){
			int size_i = graphSrc[i + 1] - graphSrc[i];
                	int size_j = graphSrc[j + 1] - graphSrc[j];
			if(size_i>size_j)
				A[j][i]=getWeight(j,i);//given j,compute i
			else
				A[i][j]=getWeight(i,j);//given i,compute j
		}
	//node n denote the '#'
	for(int i=0;i<maxVertexId;i++){
		int deg=graphSrc[i+1]-graphSrc[i];
		if(deg!=0)
			deg--;//to compute n elements ,we only need n-1 additions
		A[maxVertexId][i]=deg;
		}
	#ifdef DEBUG
	for(int i=0;i<maxVertexId+1;i++){
		for(int j=0;j<maxVertexId+1;j++)
			printf("%d	",A[i][j]);
		printf("\n");
	}
	printf("weight matrix is finished\n");
	#endif
	//generate mst store in mstSrc, mstDst, mstWeight
		
	int len=maxVertexId+1;
	int srcId=maxVertexId;
	int minId=-1;
	int minWeight=100000;
	bool* isvisited = new bool[len];
	memset(isvisited,0,sizeof(bool)*len);
	int* lowcost=new int[len];
	int* mst=new int[len];
	// memset(mst, srcId, sizeof(int)*(len));
	for(int i=0;i<len;i++){
		mst[i]=srcId;
		lowcost[i]=A[srcId][i];//start from node '#',whose ID is maxVertexId
	}
	isvisited[srcId]=true;
	lowcost[srcId]=0;
	for(int i=1;i<len;i++){
		//choose one edge to insert
		minWeight=100000;
		for(int j=0;j<len;j++){
			if(minWeight>lowcost[j]&&isvisited[j]==false){
				minId=j;
				minWeight=lowcost[j];
			}
		}
		//update isvisited, lowcost,mst
		isvisited[minId]=true;
		mstSrc[i-1]=mst[minId];
		mstDst[i-1]=minId;
		mstWeight[i-1]=A[mstSrc[i-1]][mstDst[i-1]];
		for(int j=0;j<len;j++){
			if(lowcost[j]>A[minId][j]){
				lowcost[j]=A[minId][j];
				mst[j]=minId;
			}
		}
	}
	//output mst tree
	#ifdef dEBUG
	printf("output mst\n");
	for(int i=0;i<maxVertexId;i++)
		printf("edge %d :%d %d %d\n",i,mstSrc[i],mstDst[i],mstWeight[i]);
        #endif
	delete []lowcost;
	delete []mst;
	for(int i=0;i<maxVertexId+1;i++)
		delete [] A[i];
	delete [] A;

}
int OIPSimRank::getWeight(int a,int b){//degree[a]<degree[b]
	int bWeight=graphSrc[b+1]-graphSrc[b]-1;
	int count=0;//count a intercect with b
	set<int> myset;
	//set<int>::iterator it;
	int size_a=graphSrc[a+1]-graphSrc[a];
	for(int i=graphSrc[a];i<graphSrc[a+1];i++)
		myset.insert(graphDst[i]);
	for(int i=graphSrc[b];i<graphSrc[b+1];i++)
		if(myset.find(graphDst[i])!=myset.end())
			count++;
	int x=graphSrc[b+1]-graphSrc[b]+graphSrc[a+1]-graphSrc[a]-count-count;
	if(x>bWeight)
		x=bWeight;
	return x;
	

}
void OIPSimRank::OP(int u,int iter,double* pSum){
	double* outP = new double[maxVertexId];//outP[i] denotes outP[u][i]
	memset(outP,0,sizeof(double)*maxVertexId);
	for(int i=0;i<maxVertexId;i++){// for each edge in mst
		int w=mstSrc[i];
                int z=mstDst[i];
                int scanCost= graphSrc[w+1]-graphSrc[w] + graphSrc[z+1]- graphSrc[z];
                int addCost = graphSrc[z+1]-graphSrc[z] - 1 - mstWeight[i];
        //      double factor = 0.5;	
		if(mstSrc[i]==maxVertexId || scanCost * 0.5 > addCost){// edge directly from node '#'
			int node=mstDst[i];
			for(int k=graphSrc[node];k<graphSrc[node+1];k++)
				outP[node]+=pSum[graphDst[k]];
			if(node==u)
				srvalue[1-(iter&1)][u][node]=1;
			else	if(graphSrc[node]==graphSrc[node+1] || graphSrc[u]== graphSrc[u+1]){
					srvalue[1-(iter&1)][u][node]=0;
				}
				else{
					srvalue[1-(1&iter)][u][node] = outP[node]* decayFactor / (graphSrc[u+1]-graphSrc[u])/(graphSrc[node+1]-graphSrc[node]);
				}
		}
		else {//given w, compute z
		//	int w=mstSrc[i];	
		//	int z=mstDst[i];
			outP[z]=outP[w];	
			//dst[i] of every vertex is sorted.
	        	int ws=graphSrc[w],we=graphSrc[w+1];
			int zs=graphSrc[z],ze=graphSrc[z+1];
			while(true){
				if(ws==we){
                                        for(int x=zs;x<ze;x++)
                                                outP[z]+=pSum[graphDst[x]];
                                        break;
                                        }
                                if(zs==ze){
                                        for(int x=ws;x<we;x++)
                                                outP[z]-=pSum[graphDst[x]];
                                        break;
                                        }
				if(graphDst[ws]<graphDst[zs]){
					outP[z]-=pSum[graphDst[ws]];
					ws++;
					}
				else if(graphDst[ws]>graphDst[zs]){
					outP[z]+=pSum[graphDst[zs]];
					zs++;
					}
				     else {
						zs++;ws++;
					  }
			}	
		
			if(u==z)
				srvalue[1-(iter&1)][u][z]=1;
			else   if(graphSrc[u]==graphSrc[u+1] || graphSrc[z]== graphSrc[z+1])
                        		srvalue[1-(iter&1)][u][z]=0;
                        	else
                        		srvalue[1-(1&iter)][u][z] = outP[z]* decayFactor / (graphSrc[u+1]-graphSrc[u])/(graphSrc[z+1]-graphSrc[z]);
	
		}
	
	}
	delete []outP;
}
