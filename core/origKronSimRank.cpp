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
#include "origKronSimRank.h"
void origKronSimRank::run(int qv,int k){
	if(isInit==false){//generate U,A,Vr
		isInit=true;
		firstRun=false;
		initialize();
	}
	else{//read Ku,Kv,A,Vr from file
	//...to write
		if(firstRun==true){
			A.load(Apath);
			V_r.load(V_rpath);		
			Ku.load(Kupath);
			Kv.load(Kvpath);
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
/*
void KronSimRank::initialize(){
	mat W_t(maxVertexId,maxVertexId,fill::zeros);
	for(int i=0;i<maxVertexId;i++){
		int e= origGraphSrc[i+1], s=origGraphSrc[i];
		for(int j=s;j<e;j++){
			W_t(i,origGraphDst[j])=1.0;	
		}
	}
	printf("W_t is read\n");
	double sumRow;
	double Degr=0;
	for(int i=0;i<maxVertexId;i++){
		sumRow=0;
		for(int j=0;j<maxVertexId;j++)
			sumRow+=W_t(i,j);
		if(sumRow==0.0)
			printf("node %d has no ingoing edges\n",i);
		else{
			for(int j=0;j<maxVertexId;j++)
				W_t(i,j)=W_t(i,j)/sumRow;
		}
		Degr+=sumRow;
	}
	printf("total degree is %lf\n",Degr);

	//start SVD for W_t
	Col<double> s;
	Mat<double> V_t;
	svd(U,s,V_t,W_t);//finish at 15min
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
	A.save(Apath);
	V_r.save(V_rpath);
	U.save(Upath);
	}
}
*/

void origKronSimRank::initialize(){
	mat W_t(maxVertexId,maxVertexId,fill::zeros);
        for(int i=0;i<maxVertexId;i++){
                int e= origGraphSrc[i+1], s=origGraphSrc[i];
                for(int j=s;j<e;j++){
                        W_t(i,origGraphDst[j])=1.0;     
                }
        }
	#ifdef DEBUG
        printf("W_t is read\n");
	#endif
        double sumRow;
        double Degr=0;
        for(int i=0;i<maxVertexId;i++){
                sumRow=0;
                for(int j=0;j<maxVertexId;j++)
                        sumRow+=W_t(i,j);
                if(sumRow==0.0){
			#ifdef DEBUG
                        printf("node %d has no ingoing edges\n",i);
			#endif
		}
                else{
                        for(int j=0;j<maxVertexId;j++)
                                W_t(i,j)=W_t(i,j)/sumRow;
                }
                Degr+=sumRow;
        }
	#ifdef DEBUG
        printf("total degree is %lf\n",Degr);
	#endif
	//start svd
		Mat<double> U;
		Col<double> s;
		Mat<double> V_t;
		svd(U,s,V_t,W_t);
//		s.print();
		//diagmat(s).print();
		mat V=V_t.t();

		U=U.submat(0,0,maxVertexId-1,Rank-1);
		V=V.submat(0,0,Rank-1,maxVertexId-1);
		s=s.submat(0,0,Rank-1,0);
		Ku=kron(U,U);//kronecker roduct
		mat sigma=kron(s,s);//one column
		mat K_sigma=diagmat(sigma);
		Kv=kron(V,V);

		mat K_vu=Kv*Ku;

		mat I(maxVertexId,maxVertexId);
		I.eye();
		A=inv(inv(K_sigma)-decayFactor*K_vu);
		V_r=Kv*vectorise(I);
//save temp matrix
		A.save(Apath);
		V_r.save(V_rpath);
		Kv.save(Kvpath);
		Ku.save(Kupath);
//		((1-DECAY_FACTOR)*(vectorise(I)+DECAY_FACTOR*K_u*A*V_r)).print();
}


double origKronSimRank::getScore(int i,int j){//if i==j, then return 0;
		mat V_l=(Ku.row(i*maxVertexId+j))*A;
		Mat<double> temp;
		if(i ==j)
			return 0;
		else{
			temp=(1-decayFactor)*(0+decayFactor*V_l*V_r);
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
