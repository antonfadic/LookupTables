#include "stdlib.h" // for dynamic allocation
#include "math.h" // math functions
#include "stdio.h" // in out

double *exactFun(double *inVal, int nDimOut, int nDimIn){
    //this is the exact function to compare to.
    double *y=0; y = (double*) calloc(sizeof(double),nDimIn);
    double x[nDimIn];
    int i=0;
    for (i=0;i<nDimIn;i++){x[i] = *(inVal+i);}

    *y = pow(x[0],5); //first output
    return y;
}

double *exactDeriv(double *inVal, int outDim, int nDimIn){
    //This returns the partial derivatives
    int i;
    double x[nDimIn]; for(i=0;i<nDimIn;i++){x[i]=*(inVal+i);}

    double *out; out = (double*) calloc(sizeof(double),nDimIn);

    *(out+0) = 4*pow(x[0],3); //partial derivative with respect to x[0]
    *(out+1) = 0; //partial derivative with respect to x[1]
    //*(out+2) = 0; //partial derivative with respect to x[2]
    //*(out+3) = 0; //partial derivative with respect to x[3]

    return out;
}
