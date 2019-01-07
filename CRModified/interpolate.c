#include "hInterp.h" // include exact function
#include "stdlib.h" // for dynamic allocation
#include "math.h" // math functions

double *interpolate(double *x, double *tabConfig, double *ptrTableFirst){
    //this function interpolates values one by one.
    int i,k;
    int nDimIn = *(tabConfig); // 0th parameter of table config
    int nDimOut = *(tabConfig+1); //number of dimensions in. 1st parameter
    double *nBreaks; nBreaks = tabConfig+2; //number of breaks: third parameter of the pointer ptrTableFirst
    double *min; min = (tabConfig + 2 + nDimIn);
    double *max; max = (tabConfig + 2 + 2*nDimIn);

    int nVals = getNumVals(tabConfig+2, nDimIn);

    double *h = 0; h = (double*) malloc(sizeof(double)*nDimIn); //h is the spacing in between the data.
    for(i=0;i<nDimIn;++i){ *(h+i) = (*(max+i)- *(min+i))/(*(nBreaks+i)-1); } //get the h
    double *t=0; t=(double*) malloc(sizeof(double)*nDimIn); //position in table
    int *fl=0; fl=(int*) malloc(sizeof(int)*nDimIn);
    double *xTemp=0; xTemp=(double*) malloc(sizeof(double)*nDimIn);
    for(i=0;i<nDimIn;++i){*(xTemp+i)= *(x+i);}; //copy query vector

    //calculate t
    for(k=0;k<nDimIn;k++){
        *(fl+k)= (int) floor((*(xTemp+k) - *(min+k)) / (*(h+k))) -1; //actually it is fll
        *(t+k)= (*(xTemp+k)- *(min+k))/(*(h+k)) - (double) *(fl+k)-1;
       // printf("dim %i, x is %f, h %f, fl %i and t is %f \n", k+1, *(xTemp+k) , *(h+k),*(fl+k), *(t+k));
    }

    int *pos=0;
    pos = getPosition(fl, nBreaks, nDimIn); //get the position if it lies on the borders for interpolation order
    //checkInBoundaries(fl,t,nDimIn,nBreaks);       //adjust position if it lies outside borders. Safety measure to avoid extrapolation.

    //get the 4D hypercube.
    double *answer; answer = getAnswer(nDimIn, ptrTableFirst, fl, nDimOut,nBreaks, nVals);

    //Combined Method
    //answer = CombiMultiDim(nDimIn,nDimOut,answer,t,pos); //check left. Accurate but fails - check right. No fail, not accurate.
    //Pure CR
    //answer = CRMultiDim(nDimIn,nDimOut,answer,t,pos); // Min 5 breaks.
    //PolyMultiDim
    answer = PolyMultiDim(nDimIn,nDimOut,answer,t,pos); //check right. Does not converge.

    free(xTemp); free(pos); free(fl); free(t); free(h);

    return answer;
}


