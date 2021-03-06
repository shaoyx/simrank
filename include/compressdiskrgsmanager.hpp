#ifndef __COMPRESSDISKRGSMANAGER_H__
#define __COMPRESSDISKRGSMANAGER_H__

#include "config.h"
#include "gsinterface.h"
#include "compressdiskbasedrsg.hpp"

class CompressDiskRGSManager : public GSInterface{
public:
	CompressDiskRGSManager(int sn, int mvid, char* prefix) :
		sampleGraphNum(sn), maxVertexId(mvid), sampleId(-1) {
		drsg = new CompressDiskBasedRSG(prefix);
	}

	~CompressDiskRGSManager(){
		delete drsg;
	}

	void initialize(){
		drsg->initialize(maxVertexId);
	}

	void insertEdge(int sid, int src, int dst){
//        printf("\t%d %d %d\n", sid, src, dst);
		if(sampleId == -1) sampleId = sid;
		if(sid != sampleId){
			drsg->setNewSampleId(sampleId);
			drsg->preprocess();
			drsg->save();
			drsg->clear();
			sampleId = sid;
		}
		/*revserse the edge (src, dst) here. */
		drsg->addEdge(dst, src);
	}

	void save(int sid){
		drsg->setNewSampleId(sampleId);
		drsg->preprocess();
		drsg->save();
	}

	void computeSimrank(int sid, vector<SimRankValue>& sim, map<int, vector<int>* >& timestamp, int maxSteps, double df, int qv){
		double buildCost = 0.0;
        double precost = 0.0;
        double sortcost = 0.0;
		Time timer;
		timer.start();
		/* 1. reverse graph */
		map<int, vector<int>* >::iterator iter;
		for(iter = timestamp.begin(); iter != timestamp.end(); ++iter){
			sort((*(iter->second)).begin(), (*(iter->second)).end());
		}
        timer.stop();
        sortcost = timer.getElapsedTime();

		/* load sample graph here. */
		drsg->setNewSampleId(sid);
		drsg->read();
        timer.stop();
        precost = timer.getElapsedTime();

		/* 2. compute simrank. */
		int comp = 0;
		int traverse_cnt = 0;
		/* enumerate meeting points here. */
		for(iter = timestamp.begin(); iter != timestamp.end(); ++iter){
			vector<int>* tsv = iter->second;
			int vid = iter->first; /* one meeting points. */

			int idx = 0;
			int stepLim;
//			vector<int> cand[2];
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
//					sort(cand[step&1].begin(), cand[step&1].end());
					//int size = cand[step & 1].size();
//					cand[(step + 1)&1].clear();
//					int curNbrVid = 0; /* indicate for scanning the sample graph */

//					for(int zx = 0; zx < size; ++zx){
					while(cand[step&1].empty() == false){
						cur = cand[step & 1].front(); //cand[step & 1][zx];
						cand[step & 1].pop();
						traverse_cnt++;
                        
         //               printf("step=%d: cur=%d\n", step, cur);

						if(cur != qv && step == stepLim){
							comp += cnt;
							/* update SimRank */
							sim[cur].setVid(cur);
							/* non-first meeting guarantee. */
							sim[cur].incValue(pow(df, step)*cnt);
						}

						/* enumerate edge here! */
						drsg->expand(cur, cand[(step + 1) & 1]);//, curNbrVid);
//                        int nLen = reversedSampleGraph[sid][cur].size();
//						for(int i = 0; i < nLen; ++i){
//							cand[(step + 1) & 1].push(reversedSampleGraph[sid][cur][i]);
//						}
					}
					step++;
				}while(step <= stepLim);
                //printf("stepLim=%d step=%d cnt=%d\n", stepLim, step, cnt);
			}
		}
		timer.stop();
		buildCost = timer.getElapsedTime();
        printf("sid=%d sim_comp=%d real_steps=%d comp_cost=%.5lf pre_cost=%.5lf sort_cost=%.5lf\n", sid, comp, traverse_cnt, buildCost, precost, sortcost);

	}

private:
	int sampleGraphNum;
	int maxVertexId;
	int sampleId;
	CompressDiskBasedRSG* drsg;
};

#endif
