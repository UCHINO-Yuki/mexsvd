#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))
void sgesvdx(mxChar *JOBU, mxChar *JOBVT, mxChar *RANGE, mwSignedIndex *m, mwSignedIndex *n, float *A, mwSignedIndex *LDA, float *VL, float *VU, mwSignedIndex *IL, mwSignedIndex *IU, mwSignedIndex *NS, float *S, float *U, mwSignedIndex *LDU, float *V, mwSignedIndex *LDV, float *work, mwSignedIndex *lwork, mwSignedIndex *iwork, mwSignedIndex *info);
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	mwSignedIndex i, m, n, MAT, lwork, info, IL, IU, NS;
    mwSignedIndex *iwork;
    mxChar *JOBU, *JOBVT, RANGE;
	float *A, *B, *U, *S, *V, *Stmp, *work, VL, VU;
    float dammy, *res;

    // Arguments
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	A = (float *)mxGetPr(prhs[0]);
    JOBU = mxGetChars(prhs[1]);
    JOBVT = mxGetChars(prhs[2]);
    RANGE = 'A';
    MAT = (mwSignedIndex)mxGetScalar(prhs[3]);
    iwork = MALLOC(mwSignedIndex, n*12);
    Stmp = MALLOC(float, n);
    B = MALLOC(float, m*n);

    // Set Arguments
    if (nlhs == 1) {
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
            plhs[1] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
            S = (float *)mxGetPr(plhs[1]);
        } else {
            plhs[1] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
            S = (float *)mxGetPr(plhs[1]);
        }
        plhs[0] = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
        U = (float *)mxGetPr(plhs[0]);
        if (*JOBVT != 'N') {
            plhs[2] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
            V = (float *)mxGetPr(plhs[2]);
        } else {
            V = MALLOC(float, n*n);
        }
    }

    for (i=0; i<m*n; i++) B[i] = A[i];

    // Get the length
    lwork = -1;
    sgesvdx(JOBU, JOBVT, &RANGE, &m, &n, B, &m, &VL, &VU, &IL, &IU, &NS, Stmp, U, &m, V, &n, &dammy, &lwork, iwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy;
    work = MALLOC(float, lwork);

    // Run SVD
    sgesvdx(JOBU, JOBVT, &RANGE, &m, &n, B, &m, &VL, &VU, &IL, &IU, &NS, Stmp, U, &m, V, &n, work, &lwork, iwork, &info);
    
    // Set the return values
    if (MAT == 1) {
        for (i=0; i<n; i++) S[i*n+i] = Stmp[i];
    } else {
        for (i=0; i<n; i++) S[i] = Stmp[i];
    }

    // Finalize
    free(Stmp);
    free(work);
    free(iwork);
    free(B);
    if (*JOBU == 'N') free(U);
    if (*JOBVT == 'N') {
        free(V);
        if (nlhs == 3) {
            plhs[2] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
            res = (float *)mxGetPr(plhs[2]);
            *res = 0;
        }
    }
}
