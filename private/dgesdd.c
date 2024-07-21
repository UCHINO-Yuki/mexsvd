#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))
void dgesdd(mxChar *JOBZ, mwSignedIndex *m, mwSignedIndex *n, double *A, mwSignedIndex *LDA, double *S, double *U, mwSignedIndex *LDU, double *V, mwSignedIndex *LDV, double *work, mwSignedIndex *lwork, mwSignedIndex *iwork, mwSignedIndex *info);
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	mwSignedIndex i, m, n, MAT, ecoS, lwork, info;
    mwSignedIndex *iwork;
    mxChar *JOBZ;
	double *A, *B, *U, *S, *V, *Stmp, *work;
    double dammy[1];

    // Arguments
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	A = mxGetPr(prhs[0]);
    JOBZ = mxGetChars(prhs[1]);
    MAT = (mwSignedIndex)mxGetScalar(prhs[2]);
    ecoS = (mwSignedIndex)mxGetScalar(prhs[3]);
    Stmp = MALLOC(double, n);
    B = MALLOC(double, m*n);
    iwork = MALLOC(mwSignedIndex, 8*n);

    // Set Arguments
    if (nlhs == 1) {
        if (*JOBZ == 'S') {
            if (MAT == 0) {
                plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
                S = mxGetPr(plhs[0]);
            } else {
                plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
                S = mxGetPr(plhs[0]);
            }
            U = MALLOC(double, m*n);
            V = MALLOC(double, n*n);
        } else {
            if (MAT == 0) {
                plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
                S = mxGetPr(plhs[0]);
            } else {
                if (ecoS == 0) plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
                else plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
                S = mxGetPr(plhs[0]);
            }
            U = MALLOC(double, m*m);
            V = MALLOC(double, n*n);
        }
    } else {
        if (*JOBZ == 'S') {
            if (MAT == 0) {
                plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
                S = mxGetPr(plhs[1]);
            } else {
                plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
                S = mxGetPr(plhs[1]);
            }
            plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
            U = mxGetPr(plhs[0]);
            plhs[2] = mxCreateDoubleMatrix(n, n, mxREAL);
            V = mxGetPr(plhs[2]);
        } else {
            if (MAT == 0) {
                plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
                S = mxGetPr(plhs[1]);
            } else {
                plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
                S = mxGetPr(plhs[1]);
            }
            plhs[0] = mxCreateDoubleMatrix(m, m, mxREAL);
            U = mxGetPr(plhs[0]);
            plhs[2] = mxCreateDoubleMatrix(n, n, mxREAL);
            V = mxGetPr(plhs[2]);
        }
    }

    for (i=0; i<m*n; i++) B[i] = A[i];

    // Get the length
    lwork = -1;
    dgesdd(JOBZ, &m, &n, B, &m, Stmp, U, &m, V, &n, dammy, &lwork, iwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);

    // Run SVD
    dgesdd(JOBZ, &m, &n, B, &m, Stmp, U, &m, V, &n, work, &lwork, iwork, &info);
    
    // Set the return values
    if (MAT == 1) {
        if (ecoS == 1) {
            for (i=0; i<n; i++) S[i*n+i] = Stmp[i];
        } else {
            for (i=0; i<n; i++) S[i*m+i] = Stmp[i];
        }
    } else {
        for (i=0; i<n; i++) S[i] = Stmp[i];
    }

    // Finalize
    free(Stmp);
    free(work);
    free(iwork);
    free(B);
    if (nlhs == 1) {
        free(U);
        free(V);
    }
}
