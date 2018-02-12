/*************************************************************************
    > File Name: KroSimRank.cpp
    > Author: ma6174
    > Mail: ma6174@163.com 
    > Created Time: Tue 10 Mar 2015 05:58:41 PM CST
 ************************************************************************/
/*
 * input is the normalized adjacency matrix W_t
 * output is the similarity matrix
 * the matrix V in this method is V_T traditionally
 */
//#include "/home/zippo/matrixCode/include/mybase.h"
//using namespace arma;
#include "KronSimRank.h"
void KronSimRank::run(int qv,int k){
	if(isInit==false){//generate U,A,Vr
		isInit=true;
		firstRun=false;
		initialize();
	}
	else{//read U,A,Vr from file
	//...to write
		if(firstRun==true){
			A.load(Apath);
			V_r.load(V_rpath);		
			U.load(Upath);	
		}	
	}
	double score=0.0;
	vector<SimRankValue> res;
	res.resize(maxVertexId);
	for(int i=0;i<maxVertexId;i++){
		if(qv == i){
			score = 0;//the topk cannot be itself
		}
		else {
			score=getScore(qv,i);
		}
		res[i].setVid(i);
		res[i].setValue(score);
	}	
	save(res,k);
	//compute simrank value for(S(qv,other))
	
}
void KronSimRank::initialize(){
//	printf("into initialize maxVertexid=%d\n",maxVertexId);
//	for(int i=0;i<=maxVertexId;i++)
//		printf("origGraphSrc[%d]:%d\n",i,origGraphSrc[i]);
	mat W_t(maxVertexId,maxVertexId,fill::zeros);
	for(int i=0;i<maxVertexId;i++){
//		printf("dealing with %dth vertex\n",i);
		int e= origGraphSrc[i+1], s=origGraphSrc[i];
//		printf("s is %d, e is %d\n",s,e);
		for(int j=s;j<e;j++){
			W_t(i,origGraphDst[j])=1.0;	
		}
	}
	W_t = W_t.t();
	#ifdef DEBUG
	printf("W_t is read\n");
	#endif
	double sumRow;
	double Degr=0;
	//***********************************
	vector<double> d;
	d.resize(maxVertexId);
	//***********************************
	for(int i=0;i<maxVertexId;i++){
		sumRow=0;
		for(int j=0;j<maxVertexId;j++)
			sumRow+=W_t(i,j);
		if(sumRow < 0.00000000001){
			#ifdef DEBUG
			printf("node %d has no ingoing edges\n",i);
			#endif
		//	d[i] = 1 - decayFactor;
			d[i] = 1;
		}
		else{
			for(int j=0;j<maxVertexId;j++){
				W_t(i,j)=W_t(i,j)/sumRow;
				d[i] = (sumRow - 1.0) / sumRow;
		//		d[i] = 1;
		//		d[i] = 1 - decayFactor;
				
			}
		}
		Degr+=sumRow;
	}
	#ifdef DEBUG
		printf("total degree is %lf\n",Degr);
	#endif
	//start SVD for W_t
	Col<double> s;
	Mat<double> V_t;
	svd(U,s,V_t,W_t);//finish at 15min
	#ifdef DEBUG
		printf("SVD finsihed\n");
	#endif
	mat V=V_t.t();
	U=U.submat(0,0,maxVertexId-1,Rank-1);
	#ifdef DEBUG
		printf("sub U got\n");
	#endif
	V=V.submat(0,0,Rank-1,maxVertexId-1);
	mat vu=V*U;
	mat K_vu=kron(vu,vu);
	#ifdef DEBUG
		printf("K_vu got\n");
	#endif
	s=s.submat(0,0,Rank-1,0);
	mat sigma=kron(s,s);//one column
	mat sigma_1=1./sigma;
	mat K_sigma_1=diagmat(sigma_1);
	#ifdef DEBUG
		printf("K_sigma_-1 got\n");
	#endif
	mat I(maxVertexId,maxVertexId);
	I.eye();
	A=inv(K_sigma_1-decayFactor*K_vu);
	#ifdef DEBUG
	printf("A got\n");//finish at 48min	
	#endif
//	V_r=kron(V,V)*vectorise(I);// this needs too much memory
				
	//need to be rewrite
	V_r = d[0] * kron(V.col(0),V.col(0));
	for(int i=1;i<maxVertexId;i++){
		V_r = V_r + d[i] *kron(V.col(i),V.col(i));
	}
//	char Upath[125];
//	char V_rpath[125];
//	char Apath[125];
//	sprintf(Upath,"dataset/%s/OptKron/%s.U",graphName,graphName);
//	sprintf(V_rpath,"dataset/%s/OptKron/%s.V_r",graphName,graphName);
//	sprintf(Apath,"dataset/%s/OptKron/%s.A",graphName,graphName);
		A.save(Apath);
		V_r.save(V_rpath);
		U.save(Upath);
}
/*
void KronSimRank::preprocess(Mat<double> &W_t){
		Mat<double> U;
		Col<double> s;
		Mat<double> V_t;
		svd(U,s,V_t,W_t);
//		s.print();
		//diagmat(s).print();
		mat V=V_t.t();

		U=U.submat(0,0,NODE-1,SVD_RANK-1);
		V=V.submat(0,0,SVD_RANK-1,NODE-1);
		s=s.submat(0,0,SVD_RANK-1,0);
		K_u=kron(U,U);//kronecker roduct
		mat sigma=kron(s,s);//one column
		mat K_sigma=diagmat(sigma);
		K_v=kron(V,V);

		mat K_vu=K_v*K_u;

		mat I(NODE,NODE);
		I.eye();
		A=inv(inv(K_sigma)-DECAY_FACTOR*K_vu);
		V_r=K_v*vectorise(I);
//		((1-DECAY_FACTOR)*(vectorise(I)+DECAY_FACTOR*K_u*A*V_r)).print();
}
*/

double KronSimRank::getScore(int i,int j){
		mat V_l=kron(U.row(i),U.row(j))*A;
		Mat<double> temp;
		if(i == j)
			return 0;
		else{
	//		temp=(1-decayFactor)*(0+decayFactor*V_l*V_r);
			temp = decayFactor * V_l * V_r;
			return temp(0,0);
		}
}

/*
void KronSimRank::topK(int q,int topk){
	double theta=0;
	for(int i=0;i<topk;i++){
		coll.push(make_pair(0,i));
	}
	double S_qi;
	for(int i=0;i<NODE;i++){
	if(i==q) continue;
	S_qi=getScore(q,i);
	if(S_qi>theta){
		coll.pop();
		coll.push(make_pair(S_qi,i));
		theta=coll.top().first;
	}	
	}
}
*/
