/* GKMKernel.cpp : modified GKMKernel function for python interfacing
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 * modified by 2020 Seong Kyu Han
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>       // std::cout
#include <thread>         // std::thread

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <iostream>

#include "global.h"
#include "globalvar.h"

#include "Sequence.h"
#include "CountKLmers.h"
#include "CountKLmersGeneral.h"
#include "CountKLmersH.h"
#include "CalcWmML.h"
#include "MLEstimKLmers.h"
#include "MLEstimKLmersLog.h"
#include "KLmer.h"
#include "SequenceNames.h"
#include "EstimLogRatio.h"
#include "LTree.h"
#include "LTreef.h"
#include "LTreeS.h"
#include "LList.h"
#include "CiDLPasses.h"
#include "GTree2.h"

#include "GkmKernel.h"
#include "clog.h"

using namespace std;

double calcinnerprod(int i, int j, double *c)
{
	double res = 0; 
	for(int m=0;m<=::gMAXMM;m++)
	{
		res+=::gMMProfile[i][m][j]*c[m]; 
	}
	return(res); 
}
// gives inner prodict of the pseudo-counts . nA is the number of L-mers in A and is equal to length(A)-L+1, btL is b^L
double calcinnerprod(int i, int j, double *c, double n0, double C, int nA, int nB, double btL)
{
	double res = 0; 
	for(int m=0;m<=::gMAXMM;m++)
	{
		res+=::gMMProfile[i][m][j]*c[m]; 
	}

	res = res+(nA+nB)*n0*C+btL*n0*n0; 
	return(res); 
}

void task1(int L, int j0, CiDLPasses *iDL, CLTreeS *seqsTS, int M, int nThreads)
{
	//    sprintf(globtmpstr,"started pass %d out of %d.\n",j+1,iDL->M);Printf(globtmpstr);
	for(int j=0;j<M;j++) if((j%nThreads)==j0){
		
		printf("Thread %d, started pass %d out of %d.\n",j0+1, (j-j0)/nThreads+1,1+(M-j0-1)/nThreads);
		
		int *tmpArray1 = new int[L];
		int *tmpArray2 = new int[L];
		//for(int j=0;j<iDL->M;j++){
		CLTreeS *seqsTSj= new CLTreeS();
		seqsTS->cloneReorder(seqsTSj, iDL->passOrder[j], L,L,globalConverter.b, tmpArray1, tmpArray2);
		//seqsTS->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
		//gDFSlistT[0][0] = seqsTSj; // with nonEmptyDaughterCnt
		//gDFSMMlist[0][0] = 0;
		//if(!((iDL->passTrees[j]->child0==NULL)&&(iDL->passTrees[j]->child1==NULL))) // i.e. if not empty tree
		//    seqsTSj->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL->passTrees+j, 0, globalConverter.b);
		int zero=0;
		if(!((iDL->passTrees[j]->child0==NULL)&&(iDL->passTrees[j]->child1==NULL))) // i.e. if not empty tree
			seqsTSj->DFSTiDL(&seqsTSj,1, &zero, iDL->passTrees+j, 0, globalConverter.b);
		seqsTSj->deleteTree(L, globalConverter.b, 1);
		delete seqsTSj;
		
		// print mismatch profile:
		/* for(int si=0;si<nseqs; si++){
		for(int sj=0;sj<nseqs;sj++){
		printf("\n (s%d, s%d) = ",si,sj);
		for(int dd = 0; dd<=gMAXMM; dd++){
		printf("%d ",gMMProfile[si][dd][sj]);
		}
		}
	}
		*/
		
		
		//}
		delete []tmpArray1;
		delete []tmpArray2;
		//sprintf(globtmpstr,"ended pass %d out of %d.\n",j+1,iDL->M);Printf(globtmpstr);
		}
	
	}

// Maing Kernel
void gkmKernelSuffixTree(OptsGkmKernel &opt, double **kmat, int *narr)
{
	
	int L = opt.L;
	int K = opt.K;
	int maxseqlen =	opt.maxseqlen;
	int useTgkm = opt.useTgkm;
	int maxnmm = opt.maxnmm; //auto 
	int nMAXSEQUENCES = opt.maxnumseq; 
	bool addRC = opt.addRC;
	bool usePseudocnt= opt.usePseudocnt; 
	
	char *posSeqsFN = opt.posfile;
	char *negSeqsFN = opt.negfile;
	
	int i = 1; 
	
	//char tmps[1000];
	
	char **seqname = new char *[nMAXSEQUENCES];
	
	CCalcWmML wmc(L, K, globalConverter.b);
	//double *kernel = wmc.kernelTruncated;
	if (maxnmm==-1)
	{ 
		maxnmm=L;
		if (useTgkm==1)
		{
			maxnmm = 2*(wmc.kernelTruncatedLength-1);
			if (maxnmm>L) 
			{
				maxnmm=L;
			}
		}
		if (useTgkm==2)
		{
			maxnmm = L-K; 
		}		
		if (useTgkm==3)  //wildcard kernel 
		{
			maxnmm = opt.wildcardMismatchM;
		}
		if (useTgkm==4)  //mismatch kernel 
		{
			maxnmm = 2*opt.wildcardMismatchM;
		}
		if (maxnmm>L)
		{
			maxnmm=L;
		}
	}
	double n0 = wmc.n0; 
	double *c = wmc.cTr; 
	
	n0 = c[maxnmm]/2; 
	
	if (!useTgkm)
	{
		n0 = 0; 
		//kernel = wmc.kernel; 
		c = wmc.c; // same as kernel
	}
	if (useTgkm==2)
	{
		//	n0 = 0; 
		//	kernel = wmc.kernel; 
		c = wmc.h; 
		n0 = c[maxnmm]/2;
		
	}
	if (useTgkm==3)  //wildcard kernel 
	{
		c = wmc.calcWildcardKernelWeights(L,  opt.wildcardMismatchM, globalConverter.b, opt.wildcardLambda, c);
		n0 = c[maxnmm]/2;
		
	}
	if (useTgkm==4)  //mismatch kernel 
	{
		c = wmc.calcMismatchKernelWeights(L,  opt.wildcardMismatchM, globalConverter.b, c);
		n0 = c[maxnmm]/2;
		
	}
	
	sprintf(globtmpstr,"\n maximumMismatch = %d\n", maxnmm);Printf(globtmpstr);
	for(int ii=0;ii<=maxnmm;ii++) {
		sprintf(globtmpstr,"\n c[%d] = %e",ii,c[ii] ); 	Printf(globtmpstr);
	}
	Printf("\n");
	
	int npos=0; 
	int nneg=0;
	
	gMAXMM=maxnmm; //MaxMismatch
	int UseGTree = 0; // GTree algorithm for speed
	//if(UseGTree){
	//int maxn=(1<<(2*L)) * ::Combinations(L, gMAXMM);
	// if(maxn>100000000){maxn=100000000;};
	// int maxn=100000000;
	// gGTreeLeaves=new GTreeLeafData[maxn]; // list of all the leaf nodes
	// gGTreeLeavesCnt=0; // number of all leaf nodes
	
	//}
	
	
	CLTreeS *seqsTS= new CLTreeS();
	//    GTree *seqsGTree= new GTree();
	//   GTree2 *seqsGTree2= new GTree2();
	
	
	
	int **seqsB = new int *[nMAXSEQUENCES]; 
	int **seqsBrc  = new int *[nMAXSEQUENCES]; 
	
	int *LmersCnt = new int [nMAXSEQUENCES]; 
	
	CSequence sgii(maxseqlen+3);
	CSequence *sgi = &sgii;
	
	int ntotal = 0; //number of lmers
	int nseqs=0;
	
	char**seqname2= NULL; 
	
	seqname2 = new char *[nMAXSEQUENCES];
	
	//read positive sequence file
	FILE *sfi = fopen(posSeqsFN, "r"); 
	if (sfi == NULL)
	{
		perror ("error occurred while opening a file");
		return 0;
	}
	
	char *tmpSeq=new char[maxseqlen+3];
	int  *tmpSeqB=new int[maxseqlen+3];
	
	int multiseq_allowed=1;
	
	
	while (!feof(sfi))
	{
		
		sgi->readFsa(sfi); 
		
		if(sgi->getLength()>0)
		{
			int dupl_seq_idx= -1;
			if(multiseq_allowed){
				dupl_seq_idx = find_str(seqname,nseqs, sgi->getName());
			}
			if (dupl_seq_idx>=0){
				sgi->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
				ntotal-=LmersCnt[dupl_seq_idx];
				LmersCnt[dupl_seq_idx]+= seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
				ntotal+=LmersCnt[dupl_seq_idx];
				if(addRC)
				{
					sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
					LmersCnt[dupl_seq_idx] += seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
				}
			}else{
				
				seqname2[nseqs] = new char[100];
				sprintf(seqname2[nseqs],"%s", sgi->getName());
				seqname[nseqs]=seqname2[nseqs]; 
				
				seqsB[nseqs] = new int [sgi->getLength()]; 
				sgi->getSubseqBaseId(0, sgi->getLength()-1, seqsB[nseqs]); 
				LmersCnt[nseqs] = seqsTS->addSequence(seqsB[nseqs], sgi->getLength(),L, nseqs);
				
				if(addRC)
				{
					seqsBrc[nseqs] = new int [sgi->getLength()];
					sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, seqsBrc[nseqs]); 
					LmersCnt[nseqs] = LmersCnt[nseqs] + seqsTS->addSequence(seqsBrc[nseqs], sgi->getLength(),L, nseqs);
				}
				else
				{
					seqsBrc[nseqs]=NULL; 
				}
				
				ntotal = ntotal + LmersCnt[nseqs]; 
				nseqs++;
			}
		}
	}
	fclose(sfi);
	npos = nseqs;
	
	//read negative sequence file
	sfi = fopen(negSeqsFN, "r"); 
	while (!feof(sfi))
	{
		sgi->readFsa(sfi); 
		
		if(sgi->getLength()>0)
		{
			// check if segName already exists
			int dupl_seq_idx= -1;
			if(multiseq_allowed){
				dupl_seq_idx = find_str(seqname+npos,nseqs-npos, sgi->getName());
			}
			if (dupl_seq_idx>=0){
				dupl_seq_idx+=npos;
				sgi->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
				ntotal-=LmersCnt[dupl_seq_idx];
				LmersCnt[dupl_seq_idx] += seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
				ntotal+=LmersCnt[dupl_seq_idx];
				if(addRC)
				{
					sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
					LmersCnt[dupl_seq_idx] += seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
				}
			}else{
				seqname2[nseqs] = new char[100];
				sprintf(seqname2[nseqs],"%s", sgi->getName());
				seqname[nseqs]=seqname2[nseqs];
				
				seqsB[nseqs] = new int [sgi->getLength()];
				sgi->getSubseqBaseId(0, sgi->getLength()-1, seqsB[nseqs]);
				LmersCnt[nseqs] = seqsTS->addSequence(seqsB[nseqs], sgi->getLength(),L, nseqs);
				if(addRC)
				{
					seqsBrc[nseqs] = new int [sgi->getLength()];
					sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, seqsBrc[nseqs]);
					LmersCnt[nseqs] = LmersCnt[nseqs] + seqsTS->addSequence(seqsBrc[nseqs], sgi->getLength(),L, nseqs);
					
				}
				else
				{
					seqsBrc[nseqs]=NULL; 
				}
				
				ntotal = ntotal + LmersCnt[nseqs]; 
				nseqs++; 
				
			}
		}
	}
	fclose(sfi);
	
	delete []tmpSeq;
	delete []tmpSeqB;
	
	nneg = nseqs - npos;
	
	// global vars init: 
	gLM1=L-1;
	gMAXMM=maxnmm; //MaxMismatch
	gMMProfile=new aint **[nseqs];
	for(int i=0;i<nseqs;i++)
	{
		gMMProfile[i] = new aint*[gMAXMM+1];
		for (int j=0;j<=gMAXMM;j++)
		{
			gMMProfile[i][j]=new aint[nseqs];
			for(int k=0;k<nseqs;k++)
			{
				gMMProfile[i][j][k]=0;
			}
		}
	}
	
	int *nodesAtDepthCnt = new int[L];
	for(int i=0;i<L; i++){
		nodesAtDepthCnt[i]=0;
	}
	
	int uniqueLmerCnt = seqsTS->leavesCount(0,L, globalConverter.b, nodesAtDepthCnt);
	sprintf(globtmpstr,"\n npos %d \n nneg %d \n  ntotal %d \n nunique %d\n",npos,nneg,ntotal,uniqueLmerCnt);Printf(globtmpstr);
	int minL2 = L; if (minL2<2) minL2 = 2; 
#ifdef USE_GLOBAL
	for(int i=0;i<=minL2;i++)
	{
		//gDFSlist[i] = new LTreeSnodeData *[uniqueLmerCnt];
		//	gDFSlistT[i] = new CLTreeSptr *[uniqueLmerCnt];  // without nonEmptyDaughterCnt
		gDFSlistT[i] = new CLTreeS *[uniqueLmerCnt];  // with nonEmptyDaughterCnt
		gDFSMMlist[i] = new int[uniqueLmerCnt]; 
		gDFSMMtree[i] = new CbinMMtree *[uniqueLmerCnt];
	}
	//int *curmmcnt = gDFSMMlist[0];
	
	//gDFSlistT[0][0] = seqsTS->daughter; // without nonEmptyDaughterCnt
	gDFSlistT[0][0] = seqsTS; // with nonEmptyDaughterCnt
	gDFSMMlist[0][0] = 0; 
#endif
	//   int UseGTree = 1;
	if (!UseGTree){
		// if no IDL bound
		/*
		seqsTS->DFST(gDFSlistT[0],1, gDFSMMlist[0], 0, globalConverter.b);
		
		for(int si=0;si<nseqs; si++){
		for(int sj=0;sj<nseqs;sj++){
		printf("\n (s%d, s%d) = ",si,sj);
		for(int dd = 0; dd<=gMAXMM; dd++){
		printf("%d ",gMMProfile[si][dd][sj]);
		}
		}
		}
		*/
		// else if IDL bound then
		for(int i=0;i<L; i++){
			//     sprintf(globtmpstr,"d%d , %d\n", i, nodesAtDepthCnt[i]);Printf(globtmpstr);
		}
		
		CiDLPasses iDL;
		//iDL.newIDLPasses(L, gMAXMM);
		double p=1.0/::globalConverter.b;
		
		//iDL.initPassOrderAll(L, gMAXMM);
		//iDL.newGreedyIDLPasses(L,iDL.M,  gMAXMM, nodesAtDepthCnt, p);
		
		iDL.newGreedyIDLPasses(L,2*L,  gMAXMM, nodesAtDepthCnt, p);
		
		//iDL.newGreedy2IDLPasses(L,2*L,  gMAXMM, nodesAtDepthCnt, p);
		//iDL.newGreedy2IDLPasses(L,L,  gMAXMM, nodesAtDepthCnt, p);
		
		//iDL.newPassOrderDesignCover( L, gMAXMM, 3);// generates M passes, that gaurantee that the first k places are matches (in all the trees)
		//iDL.newGreedyIDLPasses(L,iDL.M,  gMAXMM, nodesAtDepthCnt, p);
		
		
		/*    int *tmpArray1 = new int[L];
		int *tmpArray2 = new int[L];
		for(int j=0;j<iDL.M;j++){
		sprintf(globtmpstr,"pass %d out of %d.\n",j+1,iDL.M);Printf(globtmpstr);
		CLTreeS *seqsTSj= new CLTreeS();
		seqsTS->cloneReorder(seqsTSj, iDL.passOrder[j], L,L,globalConverter.b, tmpArray1, tmpArray2);
		//seqsTS->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
		gDFSlistT[0][0] = seqsTSj; // with nonEmptyDaughterCnt
		gDFSMMlist[0][0] = 0;
		if(!((iDL.passTrees[j]->child0==NULL)&&(iDL.passTrees[j]->child1==NULL))) // i.e. if not empty tree
		seqsTSj->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
		seqsTSj->deleteTree(L, globalConverter.b, 1);
		delete seqsTSj;
		
		// print mismatch profile:
		/ * for(int si=0;si<nseqs; si++){
		for(int sj=0;sj<nseqs;sj++){
		printf("\n (s%d, s%d) = ",si,sj);
		for(int dd = 0; dd<=gMAXMM; dd++){
		printf("%d ",gMMProfile[si][dd][sj]);
		}
		}
		}
		* /
		}
		delete []tmpArray1;
		delete []tmpArray2;
		*/
		
		//myThreads[0] = std::thread(task1, L, 0, &iDL, seqsTS);
		//myThreads[1] = std::thread(task1, L, 1, &iDL, seqsTS);
		//myThreads[0].join();
		//myThreads[1].join();
		
		int nThreads=std::thread::hardware_concurrency();
		if(nThreads==0){nThreads = iDL.M;}
		if(nThreads>iDL.M){nThreads =iDL.M;}
		nThreads = iDL.M/(int)(iDL.M/nThreads);
		if(opt.maxnThread<nThreads){nThreads=opt.maxnThread;}
		
		if (nThreads<=1){
			task1( L, 0, &iDL, seqsTS, iDL.M, 1);
		}else{
			
#ifndef MULTI_THREAD_SAFE
			printf("Warning -- MULTI_THREAD_SAFE is not enabled (see src/global.h). Some values may be approximated.\n");
#endif
			
			std::thread *myThreads = new std::thread[nThreads];
			int j;
			
			for(j=0;j<nThreads;j++){
				myThreads[j] = std::thread(task1, L, j, &iDL, seqsTS, iDL.M, nThreads);
				// myThreads[j].join();
			}
			for(j=0;j<nThreads;j++){
				myThreads[j].join();
			}
			delete []myThreads;
		}
		
	}else{
		gGTreeLeaves2=new GTreeLeafData2[uniqueLmerCnt * ::Combinations(L, gMAXMM)]; // list of all the leaf nodes
		gGTreeLeavesCnt=0; // number of all leaf nodes
		GTree2 *seqsGTree2= new GTree2();
		int *tmpArray1 = new int[L];
		
		seqsTS->addToGTree(seqsGTree2, L,tmpArray1, MAX_ALPHABET_SIZE, L);
		delete []tmpArray1;
		sprintf(globtmpstr," gGTreeLeavesCnt = %d \n",gGTreeLeavesCnt);Printf(globtmpstr);
		
		for(int i=0;i<gGTreeLeavesCnt;i++){
			gGTreeLeaves2[i].process();
		}
		// normalize
		
		for(int i=0;i<nseqs;i++)
		{
			
			//            gMMProfile[i][0][i]/=Combinations(L, gMAXMM);;
			
			for (int j=0;j<=gMAXMM;j++)
			{
				//              int s=::Combinations(L-j, L-gMAXMM);
				for(int k=0;k<=i;k++)
				{
					int mmijk = (gMMProfile[i][j][k]+gMMProfile[k][j][i]);
					//                mmijk=mmijk/s;
					gMMProfile[i][j][k]=gMMProfile[k][j][i]=mmijk;
				}
			}
		}
		
		for(int i=0;i<nseqs;i++)
		{
			
			//            gMMProfile[i][0][i]/=Combinations(L, gMAXMM);;
			
			for (int j=0;j<=gMAXMM;j++)
			{
				int s=::Combinations(L-j, L-gMAXMM);
				//                int s=::Combinations(gMAXMM, j);
				for(int k=0;k<nseqs;k++)
				{
					//                    gMMProfile[i][j][k]*=2;
					int z=gMMProfile[i][j][k]/s;gMMProfile[i][j][k]=z;
				}
			}
		}
		
		
		for(int i=0;i<nseqs;i++)
		{
			
			gMMProfile[i][0][i]+=LmersCnt[i];//self;
		}
	}
	
	delete []nodesAtDepthCnt;
	
	//
#ifdef USE_GLOBAL
	for(int i=0;i<=minL2;i++)
	{
		delete []gDFSlistT[i];
		delete []gDFSMMlist[i];
		delete []gDFSMMtree[i];
	}
#endif
	/*
	for(int i=0;i<nseqs;i++)
	{
	for(int k=0;k<nseqs;k++)
	{
	for (int j=0;j<=gMAXMM;j++)
	{
	printf("(%d,%d)[%d] = %d\n",i,k,j, gMMProfile[i][j][k]);
	}
	}
	}
	*/	
	
	/// calc C 
	double C =0; 
	for(int m=0;m<=L;m++)
	{
		C+=dCombinations(L,m)*pow(1.0*globalConverter.b-1,m)*wmc.kernelTruncated[m];
		//		C+=dCombinations(L,m)*pow(3.0,m)*wmc.kernelTruncated[m];
	}
	
	double btL=pow(1.0*globalConverter.b,L);
	double *norm = new double [nseqs]; 

    for(i=0;i<nseqs;i++)
    {
        if (usePseudocnt)
        {
            norm[i] = sqrt(calcinnerprod(i,i,c,n0,C,LmersCnt[i], LmersCnt[i], btL));
        }
        else
        {
            norm[i] = sqrt(calcinnerprod(i,i,c));
        }
    }
    
    for(i=0;i<nseqs;i++)
    {	
        for(int j=0;j<nseqs;j++)
        {
            if(i>j)
            {
                if (usePseudocnt)
                {
                    kmat[i][j] = norm[i]*norm[j]<1E-50)?0.0:calcinnerprod(i,j,c, n0,C,LmersCnt[i], LmersCnt[j], btL)/(norm[i]*norm[j])
                }
                else
                {
                    kmat[i][j] = norm[i]*norm[j]<1E-50)?0.0:calcinnerprod(i,j,c)/(norm[i]*norm[j])
                }			
            }
            else if (i==j) 
            {
                kmat[i][j] = 1.0
            }
        }
    }	
	
	// pass the number of sequences to python
	narr[0] = npos
	narr[1] = nneg

	delete []norm;
	delete []LmersCnt;
	seqsTS->deleteTree(L, globalConverter.b, 0);

	for(int i=0;i<nseqs;i++)
	{
		delete []seqsB[i]; 
		if (seqsBrc[i]!=NULL) delete []seqsBrc[i]; 
		for (int j=0;j<=gMAXMM;j++)
		{
			delete []gMMProfile[i][j];
		}
		delete []gMMProfile[i];
	}
	delete []gMMProfile;
	
	delete []seqname; 
	if (seqname2!=NULL)
	{
		for(i=0;i<nseqs;i++)
		{
			delete []seqname2[i]; 
		}
		delete []seqname2; 
	}
	
	return;
}