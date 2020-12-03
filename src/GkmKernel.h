/* GKMKernel.h : modified GKMKernel function for python interfacing - header file
 *
 * Copyright (C) 2020 Seong Kyu Han
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

typedef struct {
	int L;
	int K;
	int maxnmm;
	int maxseqlen;
	int maxnumseq;
	int useTgkm;
	bool addRC;
	bool usePseudocnt;
	char *posfile;
	char *negfile;
	double wildcardLambda;  // Parameter lambda for (LK2004)
	int wildcardMismatchM;  // Parameter M for wildcard or mismatch kernels (LK2004)
	int maxnThread;         // Max number of threads
	int dummyVal;
} OptsGkmKernel;

double calcinnerprod(int i, int j, double *c, double n0, double C, int nA, int nB, double btL);
void task1(int L, int j0, CiDLPasses *iDL, CLTreeS *seqsTS, int M, int nThreads);
void gkmKernelSuffixTree(OptsGkmKernel &opt, double **kmat, int *narr);