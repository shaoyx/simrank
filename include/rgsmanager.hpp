#ifndef __RGSMANAGER_H__
#define __RGSMANAGER_H__

#include "config.h"
#include "gsinterface.h"
#include "rsamplegraph.hpp"

class RGSManager : public GSInterface{
public:
	RGSManager(int sn, int mvid) :
		sampleGraphNum(sn), maxVertexId(mvid) {
		rsg = new RSampleGraph*[sn];
		for(int i = 0; i < sn; ++i){
			rsg[i] = new RSampleGraph(maxVertexId);
		}
	}

	~RGSManager(){
        for(int i = 0; i < sampleGraphNum; ++i)
            delete rsg[i];
		delete [] rsg;
	}

	void insertEdge(int sid, int src, int dst){
		/*revserse the edge (src, dst) here. */
		rsg[sid]->addEdge(dst, src);
	}

	void analysis(){
		double totalMemCost = 0;
		for(int i = 0; i < sampleGraphNum; ++i){
			rsg[i]->preprocess();
			totalMemCost += rsg[i]->getMemSize();
		}
		printf("Index MEM cost: %.5lfMB\n", totalMemCost);
	}

	void computeSimrank(int sid, vector<SimRankValue>& sim, map<int, vector<int>* >& timestamp, int maxSteps, double df, int qv){
		double buildCost = 0.0;
		Time timer;
		/* 1. reverse graph */
		map<int, vector<int>* >::iterator iter;
		for(iter = timestamp.begin(); iter != timestamp.end(); ++iter){
			sort((*(iter->second)).begin(), (*(iter->second)).end());
		}

		timer.start();
		/* 2. compute simrank. */
		int comp = 0;
		int traverse_cnt = 0;
		/* enumerate meeting points here. */
		for(iter = timestamp.begin(); iter != timestamp.end(); ++iter){
			vector<int>* tsv = iter->second;
			int vid = iter->first; /* one meeting points. */

			int idx = 0;
			int stepLim;
			queue<int> cand[2];
			int step = 0, cur, cnt;
            int tsvLen = tsv->size();
			cand[step].push(vid);
			while(idx < tsvLen){
				stepLim = (*tsv)[idx++];
				cnt = 1;
				while(idx < tsvLen && stepLim == (*tsv)[idx]) { idx++; cnt++; }
				/* traverse the tree. */
				do{
					while(cand[step & 1].empty() == false){
						cur = cand[step & 1].front();
						cand[step & 1].pop();
						traverse_cnt++;

						if(qv != cur && step == stepLim){
							comp += cnt;
							/* update SimRank */
							sim[cur].setVid(cur);
							/* non-first meeting guarantee. */
							sim[cur].incValue(pow(df, step)*cnt);
						}

						/* enumerate edge here! */
						rsg[sid]->expand(cur, cand[(step + 1) & 1]);
//                        int nLen = reversedSampleGraph[sid][cur].size();
//						for(int i = 0; i < nLen; ++i){
//							cand[(step + 1) & 1].push(reversedSampleGraph[sid][cur][i]);
//						}
					}
					step++;
				}while(step <= stepLim);
//                printf("stepLim=%d step=%d cnt=%d\n", stepLim, step, cnt);
			}
		}
		timer.stop();
		buildCost = timer.getElapsedTime();
        //printf("sid=%d sim_comp=%d real_steps=%d comp_cost=%.5lf\n", sid, comp, traverse_cnt, buildCost);

	}

private:
	int sampleGraphNum;
	int maxVertexId;
	RSampleGraph** rsg;
};

#endif
