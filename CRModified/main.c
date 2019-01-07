/**********************************************************************
Function for building, reading, and evaluating a combined and CR lookup table.
Meant to be used as UDF in Fluent.
Compiled with GCC on windows.
Author: Anton Fadic
University of Alberta
Chemical Engineering
**********************************************************************/
//#include "udf.h" //UDF macros
#include "stdio.h" // in out
#include "stdlib.h" // for dynamic allocation
#include "math.h" // math functions
#include "hInterp.h" //calls required functions for this file
#include "exactFun.h" //calls required functions for this file
#include <unistd.h>
#include <sys/time.h>

/*********************************************************************
1) Compile with flags -o -ffast-math. Works with gcc on Win 7 (not exc).
2) Examine the derivative of combined at the boundaries.
3) Check bugs on right bound.
4) Leaking memory!
**********************************************************************/

int main(){
    int nDimIn = 6;
    int nDimOut = 1;
    int nVals;
    int i,j;
    int refIndex=4; //mesh refinement index. Exponential, don't go crazy.
    //int hermFlag =1; // 1 if writes derivative grid, 0 if not.

    long start, end;
    long avTime=0;
    struct timeval timecheck;

    srand(time(NULL));   // should only be called once
    double r[nDimIn];        // r is a random number, that is meant to be (0,1)
    for(i=0;i<nDimIn;i++){r[i]  = (double) rand()/RAND_MAX;}

    printf("Output         Exact         Error         Time ms      Size kB\n");

    double *sol=0;
    double *derTable;
    int nIt=1000; //number of repetitions of the same calculation. Done for timing purposes.


    for(j=0;j<refIndex;j++){
        writeTableConfig(nDimIn, nDimOut,5*intPow(2,j),5*intPow(2,j)); //this writes the config file for the table
        double *tabConfig; tabConfig = readTableConfig(nDimIn); //pointer to the address of the first entry of table config
        nVals = getNumVals(tabConfig+2, nDimIn); //this gets the number of values of the grid.
        //printf("Table size %f MB \n",(double) nVals*nDimOut/1000/1000*sizeof(double));

        double *grid; grid = writeGrid(tabConfig); //writes the grid. Not working for more than 8 dimension. Fragmentate memory?

        //saveGrid(grid,nDimIn,nVals); //save grid to a file. Comment if necessary
        //grid = readGrid(nDimIn,nVals); //reads the grid if stored (faster calculations)

        writeTable(nVals,nDimIn,nDimOut,grid); //write the table to a file
        //getchar();
        derTable = getDerTable(tabConfig,grid); //get derivative values. Warning, it requires tons of memory!

        free(grid); //printf("grid freed successfully '\n");

        double *ptrTableFirst; ptrTableFirst=readFile(nVals,nDimOut); // only read once. Table stored in memory ready to use.

        //for(i=0;i<nDimOut;++i){printf("table %i first value is: %f \n", i+1, *(ptrTableFirst+nVals*i));} //show table stats

        int numTestVal = 1;
        double xTest[nDimIn*numTestVal];

        for(i=0;i<nDimIn;i++){ xTest[i] = 4.33;} //1.0+(10-9)*r[i]; } // values for table

        int nTrials=1; //number of repetitions of a calculation. This is done to see if we have consistency
        //if(j==0){printf("Warning: Number of iterations is set to %i. Result for 10 000 it \n",nIt);}
        while(nTrials>0){

        gettimeofday(&timecheck, NULL);
        start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

        for(i=0;i<nIt;i++){
            // receives the pointer of the stored interpolated data.
            //sol=interpolate(xTest,tabConfig,ptrTableFirst);
            sol=interpHermite(xTest,tabConfig,ptrTableFirst,derTable);
        }
        //save the elapsed time in ms
        gettimeofday(&timecheck, NULL);
        end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec/1000;

        //code display
        printf("%f     ", *sol);
        printf("%f     ", *exactFun(xTest,nDimOut,nDimIn));
        printf("%f     ", 1*(*sol-*exactFun(xTest,nDimOut,nDimIn)));
        printf("%ld     ", (end - start));
        printf("%f\n",(double) nVals*nDimOut/1000*sizeof(double));
        avTime += (end-start); //average time

        nTrials--;
        }
    }
    printf("Average time is %ld ms \n", avTime/refIndex);

    displayMethod(sol+1); //Display the method used. This reads the value to the right of the solution
    //Display GNU version
    printf("\n\n\ngcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
    displayAuthor();
    return 0;
}
