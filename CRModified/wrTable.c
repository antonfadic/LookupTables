#include "stdio.h" // in out
#include "exactFun.h" // include exact function
#include "stdlib.h" // for dynamic allocation
#include "hInterp.h" // include exact function

void writeTableConfig(int nDimIn, int nDimOut, int brX, int brY){
    FILE * fp;
    //int i;
    double *ptr = 0;

    ptr = (double*) malloc((3*nDimIn+2)*sizeof(double)); //nBreaks per dimension
    //write number of dimensions in
    *ptr=nDimIn;
    //write the number of dimensions out
    *(ptr+1) = nDimOut;

    //number of breaks. 0-nDimIn-1
    *(ptr+2+0) = (double) brX; //brX
    *(ptr+2+1) = (double) 10;
    *(ptr+2+2) = (double) 10;
    *(ptr+2+3) = (double) 10;
    *(ptr+2+4) = (double) 7;
    *(ptr+2+5) = (double) 7;
    //*(ptr+2+6) = (double) 7;
    //*(ptr+2+7) = (double) 7;
    //*(ptr+2+8) = (double) 7;

    //Minimum and Maximum. 0-nDimIn-1
    *(ptr+nDimIn+2+0) = (double) 1; *(ptr+2*nDimIn+2+0) = (double) 10; //Temp
    *(ptr+nDimIn+2+1) = (double) 1; *(ptr+2*nDimIn+2+1) = (double) 10; //yNH3
    *(ptr+nDimIn+2+2) = (double) 1; *(ptr+2*nDimIn+2+2) = (double) 10; //yNO
    *(ptr+nDimIn+2+3) = (double) 1; *(ptr+2*nDimIn+2+3) = (double) 10; //yO2
    *(ptr+nDimIn+2+4) = (double) 1; *(ptr+2*nDimIn+2+4) = (double) 10; //yO2
    *(ptr+nDimIn+2+5) = (double) 1; *(ptr+2*nDimIn+2+5) = (double) 10; //yO2
    //*(ptr+nDimIn+2+6) = (double) 1; *(ptr+2*nDimIn+2+6) = (double) 10; //yO2
    //*(ptr+nDimIn+2+7) = (double) 1; *(ptr+2*nDimIn+2+7) = (double) 10; //yO2
    //*(ptr+nDimIn+2+8) = (double) 1; *(ptr+2*nDimIn+2+8) = (double) 5; //yO2


    fp = fopen ("tableConfig.bin", "wb");
    fwrite(ptr,sizeof(*ptr),3*nDimIn+2,fp);
    fclose (fp);
    //printf("Config file written successfully \n \n");
    free(ptr);
}

double *readFile(int length, int nDimOut){
    FILE * fp;
    double *ptrFile = 0; ptrFile = (double*) malloc(nDimOut*length*sizeof(double));
    if(ptrFile==NULL){printf("Memory not allocated. Table could not be read. Closing"); exit(0);}
    fp = fopen ("Table.bin", "rb");
    // printf("Reading table file... \n");
    fread(ptrFile,sizeof(*ptrFile),nDimOut*length,fp);
    fclose (fp);
    return ptrFile;
}

void writeTable(int length, int nDimIn, int nDimOut, double *grid){
    FILE * fp;
    int i,j;
    double *ptrTable = 0;ptrTable=(double*) malloc(nDimOut*length*sizeof(double)); if(ptrTable==NULL){printf("memory not allocated. Closing...");exit(0);}
    //double *gridTemp = 0;gridTemp=(double*) malloc(sizeof(double)*nDimIn);
    //loop for writing the table. j is number of dimensions and i is entry.
    for(i=0; i<length; ++i){
        for(j=0;j<nDimOut;++j){
            // apply inverse transformation
            // gridTemp= invFunTransformIn((grid+i*nDimIn));
            // *(ptrTable+j*length+i) = *(exactFun(gridTemp, nDimOut, nDimIn)+j);
             *(ptrTable+j*length+i) = *(exactFun(grid+i*nDimIn, nDimOut, nDimIn)+j);
        }
    }
    fp = fopen ("Table.bin", "wb");
    //printf("Writing table file... \n");
    fwrite(ptrTable,sizeof(*ptrTable),nDimOut*length,fp);
    //free(gridTemp);
    free(ptrTable);
    fclose (fp);
    //printf("table written successfully \n");
}

double *readTableConfig(int nDim){
    FILE * fp;
    double *ptr2 = 0; ptr2 = (double*) malloc((3*nDim+2)*sizeof(double)); //required number of info: nDimIn, nDimOut, 3*nDim (breaks,min,max)
    if(ptr2==NULL){printf("Memory not allocated. Closing readTableConfig"); exit(0);}
    fp = fopen ("tableConfig.bin", "rb");
    //printf("Reading config file... \n");
    fread(ptr2,1,(3*nDim+2)*sizeof(*ptr2),fp);
    fclose (fp);
    /* //uncomment if needed to check config options
    int i;
    printf("Config file read successfully \n \n"); printf("Number of input dimensions is %i \n", (int) *ptr2); printf("Number of output dimensions is %i \n", (int) *(ptr2+1));
    for(i=0;i<nDim;++i){
        printf("Dimension %i, breaks = %i, Min = %f , Max = %f \n", i+1,(int) *(ptr2+i+2), *(ptr2+nDim+i+2), *(ptr2+2*nDim+i+2));
    }*/
    return ptr2;
}

double *getDerTable(double *tabConfig, double *grid){
    int i,j;
    int nDimIn = (int) *tabConfig;
    int nDimOut = *(tabConfig+1);
    int nVals = getNumVals(tabConfig+2, nDimIn);

    double *h = 0; h = (double*) malloc(sizeof(double)*nDimIn); //h is the spacing in between the data.
    double *min=0; min = (tabConfig + 2 + nDimIn);
    double *max=0; max = (tabConfig + 2 + 2*nDimIn);
    double *nBreaks=0; nBreaks = (tabConfig+2);

    for(i=0;i<nDimIn;++i){
        *(h+i) = (*(max+i)- *(min+i))/(*(nBreaks+i)-1);
    }

    double *derTable=0; derTable = (double*) calloc(sizeof(double),nVals*nDimIn);

    for(j=0;j<nDimIn;j++){
        for(i=0;i<nVals;i++){
            *(derTable +i+j*nVals) = *(exactDeriv(grid+i*nDimIn,j,nDimIn))*(*(h+j));
        }
    }
    return derTable;
}
