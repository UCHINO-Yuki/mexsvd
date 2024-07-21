#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))
void sgesdd(mxChar *JOBZ, mwSignedIndex *m, mwSignedIndex *n, float *A, mwSignedIndex *LDA, float *S, float *U, mwSignedIndex *LDU, float *V, mwSignedIndex *LDV, float *work, mwSignedIndex *lwork, mwSignedIndex *iwork, mwSignedIndex *info);
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	mwSignedIndex i, m, n, MAT, ecoS, lwork, info;
    mwSignedIndex *iwork;
    mxChar *JOBZ;
	float *A, *B, *U, *S, *V, *Stmp, *work;
    float dammy[1];

    // Arguments
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	A = (float *)mxGetPr(prhs[0]);
    JOBZ = mxGetChars(prhs[1]);
    MAT = (mwSignedIndex)mxGetScalar(prhs[2]);
    ecoS = (mwSignedIndex)mxGetScalar(prhs[3]);
    Stmp = MALLOC(float, n);
    B = MALLOC(float, m*n);
    iwork = MALLOC(mwSignedIndex, 8*n);

    // Set Arguments
    if (nlhs == 1) {
        if (*JOBZ == 'S') {
            if (MAT == 0) {
                plhs[0] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
                S = (float *)mxGetPr(plhs[0]);
            } else {
                plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
                S = (float *)mxGetPr(plhs[0]);
            }
            U = MALLOC(float, m*n);
            V = MALLOC(float, n*n);
        } else {
            if (MAT == 0) {
                plhs[0] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
                S = (float *)mxGetPr(plhs[0]);
            } else {
                if (ecoS == 0) plhs[0] = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
                else plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
                S = (float *)mxGetPr(plhs[0]);
            }
            U = MALLOC(float, m*m);
            V = MALLOC(float, n*n);
        }
    } else {
        if (*JOBZ == 'S') {
            if (MAT == 0) {
                plhs[1] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
                S = (float *)mxGetPr(plhs[1]);
            } else {
                plhs[1] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
                S = (float *)mxGetPr(plhs[1]);
            }
            plhs[0] = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
            U = (float *)mxGetPr(plhs[0]);
            plhs[2] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
            V = (float *)mxGetPr(plhs[2]);
        } else {
            if (MAT == 0) {
                plhs[1] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
                S = (float *)mxGetPr(plhs[1]);
            } else {
                plhs[1] = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
                S = (float *)mxGetPr(plhs[1]);
            }
            plhs[0] = mxCreateNumericMatrix(m, m, mxSINGLE_CLASS, mxREAL);
            U = (float *)mxGetPr(plhs[0]);
            plhs[2] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
            V = (float *)mxGetPr(plhs[2]);
        }
    }

    for (i=0; i<m*n; i++) B[i] = A[i];

    // Get the length
    lwork = -1;
    sgesdd(JOBZ, &m, &n, B, &m, Stmp, U, &m, V, &n, dammy, &lwork, iwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);

    // Run SVD
    sgesdd(JOBZ, &m, &n, B, &m, Stmp, U, &m, V, &n, work, &lwork, iwork, &info);
    
    // Set the return values
    if (MAT == 1) {
        if (*JOBZ == 'S') {
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
