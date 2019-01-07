#include "stdio.h"
#include "stdlib.h" // for dynamic allocation
#include "math.h" // math functions
#include "hInterp.h" // include exact function

void saveGrid(double *grid, int nDimIn, int nVals){
    FILE *fp;
    fp = fopen("grid.bin","wb");
    fwrite(grid,sizeof(double),nVals*nDimIn,fp);
    fclose(fp);
}

double *readGrid(int nDimIn, int nVals){
    FILE *fp;
    fp = fopen("grid.bin","rb");
    double *grid=0; grid=calloc(sizeof(double*),nVals);
    fread(grid,sizeof(double),nVals*nDimIn,fp);
    fclose(fp);
    return grid;
}

double *writeGrid(double *ptrFirstVal){
    //This function is in its final version. Please do not touch. This function fails for more than 6 dimensions
    //because of segmentation fault. Possible workaraound: fragmentate memory.
    int i,j;
    int nDimIn = *ptrFirstVal;
    int nVals = getNumVals(ptrFirstVal+2,nDimIn);

    double *h = 0; h = (double*) malloc(sizeof(double)*nDimIn); //h is the spacing in between the data.
    double *nBreaks=0; nBreaks = (ptrFirstVal+2);
    double *min=0; min = (ptrFirstVal + 2 + nDimIn);
    double *max=0; max = (ptrFirstVal + 2 + 2*nDimIn);

    //min = funTransformIn(min);
    //max = funTransformIn(max);

    for(i=0;i<nDimIn;++i){
        *(h+i) = (*(max+i)- *(min+i))/(*(nBreaks+i)-1);
    }
    int sumNBreaks=0; for(i=0;i<nDimIn;++i){
        sumNBreaks += *(nBreaks + i);
    }

    double *rValues=0; rValues=malloc(sizeof(double)*sumNBreaks);
    /*Objective here is to write the list of possible values as rValues={X1,X2...X_b1,Y1,Y2,...Y_b2}
    with following rule:
    Xi+1=Xi+h; X0= *(min+0);
    Yi+1=Yi+h; Y0= *(min+1); */
    int counter=0; j=0;
    while(j<nDimIn){
        for(i=0;i< *(nBreaks+j);++i){
            *(rValues+counter) = *(min+j) + *(h+j)*i;
            //printf("rValues are %f \n", *(rValues + counter));
            counter++;
        } j++;
    }

    /*The objective here is to write the vectors of *grid values from the rValues vector. This is the loops that fails at nDimIn>6*/
    double *grid=0; grid=(double*) malloc(sizeof(double)*nVals*nDimIn*1); if(grid==NULL) {printf("grid memory not allocated. Closing.\n"); exit(1);}
    int temp=1; j=0; int temp2=0;

    while(j<nDimIn){
        //printf("j is %i \n", j);
        for(i=0;i<nVals;++i){
            *(grid+i*nDimIn+j) = *(rValues + temp2 + (((int) floor(i/temp)) % (int) *(nBreaks+j)) ); //hardest function ever
        }
        temp2 += (int) *(nBreaks+j);
        temp *= (int) *(nBreaks+j);
        j++;
    }

    //min = invFunTransformIn(min);
    //max = invFunTransformIn(max);

    /* //this code is for visualizing the grid. Uncomment if necessary,
    for(i=0;i<nVals;++i){j=0; while(j<nDimIn){ printf("%f ", *(grid+i*nDimIn+j)); j++;} printf("\n"); } */
        free(h); free(rValues);
    return grid;
}
