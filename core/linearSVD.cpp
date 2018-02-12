/*************************************************************************
    > File Name: linearSVD.cpp
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
#include "linearSVD.h"
void KronSimRank::run(int qv,int k){
	if(isInit==false){//generate U,A,Vr
		isInit=true;
		firstRun=false;
		initialize();
	}
	else{//read U,A,Vr from file
	//...to write
		if(firstRun==true){//to avoid the case that first run & has file V. Qsum
			Qsum.load(Qsumpath);
			V.load(Vpath);		
		}	
	}
	double score=0.0;
	vector<SimRankValue> res;
	res.resize(maxVertexId);
	for(int i=0;i<maxVertexId;i++){
		score=getScore(qv,i);
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
	mat W(maxVertexId,maxVertexId,fill::zeros);
	mat D(maxvertexId,maxVertexId,fill::zeros);
	for(int i=0;i<maxVertexId;i++){
//		printf("dealing with %dth vertex\n",i);
		int e= origGraphSrc[i+1], s=origGraphSrc[i];
//		printf("s is %d, e is %d\n",s,e);
		for(int j=s;j<e;j++){
			W(i,origGraphDst[j])=1.0;	
		}
	}
	printf("W is read\n");
	double sumRow;
	double Degr=0;
	for(int i=0;i<maxVertexId;i++){
		sumRow=0;
		for(int j=0;j<maxVertexId;j++)
			sumRow+=W(i,j);
		if(sumRow==0.0)
			printf("node %d has no ingoing edges\n",i);
		else{
			for(int j=0;j<maxVertexId;j++){
				W(i,j)=W(i,j)/sumRow;
			}
			double alpha=1.0;
			D(i,i)= alpha - decayFactor/sumRow;
		}
		Degr+=sumRow;
	}
	printf("total degree is %lf\n",Degr);
//改到了这里，还有下面的没改，linearSVD.h改完了，但是没有测试。
	//start SVD for W_t
	Col<double> s;
	Mat<double> V_t;
	svd(U,s,V_t,W);//finish at 15min
	printf("SVD finsihed\n");
	mat V=V_t.t();
	U=U.submat(0,0,maxVertexId-1,Rank-1);
	printf("sub U got\n");
	V=V.submat(0,0,Rank-1,maxVertexId-1);
	mat vu=V*U;
	mat K_vu=kron(vu,vu);
	printf("K_vu got\n");
	s=s.submat(0,0,Rank-1,0);
	mat sigma=kron(s,s);//one column
	mat sigma_1=1./sigma;
	mat K_sigma_1=diagmat(sigma_1);
	printf("K_sigma_-1 got\n");
	mat I(maxVertexId,maxVertexId);
	I.eye();
	A=inv(K_sigma_1-decayFactor*K_vu);
	printf("A got\n");//finish at 48min	
//	V_r=kron(V,V)*vectorise(I);// this needs too much memory
				
	//need to be rewrite
	V_r=kron(V.col(0),V.col(0));
	for(int i=1;i<maxVertexId;i++){
	V_r=V_r+kron(V.col(i),V.col(i));
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
		if(i!=j)
			temp=(1-decayFactor)*(0+decayFactor*V_l*V_r);
		else
			temp=(1-decayFactor)*(1+decayFactor*V_l*V_r);
		double result=temp(0,0);
		return result;
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
