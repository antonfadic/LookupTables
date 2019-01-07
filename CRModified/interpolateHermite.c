#include "hInterp.h" // include exact function
#include "stdio.h" // in out
#include "stdlib.h" // for dynamic allocation
#include "math.h" // math functions

double *interpHermite(double *x, double *tabConfig, double *ptrTableFirst,  double *derTable){
    //Hermite method
    int nDimIn = (int)  *tabConfig;
    //double *nBreaks = tabConfig+2;

    int nVals = getNumVals(tabConfig+2, nDimIn);
    int *fl = get2DCoord(x,tabConfig); //OK

    double *t = getT(x,tabConfig,fl); //bug. Wrong data type.

    //Get the 2D hypercube. Bug
    double *funTable = get2DFun(fl, tabConfig, ptrTableFirst);

    //get the 2D Derivative table
    double *derAnswer = 0;
    derAnswer = get2DderTable(fl, tabConfig, derTable);


    int offset=0;
    int j=0; int i=0;


   while(j<nDimIn){
        //read the derivatives

        for(i=0;i<intPow(2,nDimIn-j-1);i++){
            *(funTable+i) = markHermite(*(derAnswer+2*i),*(funTable+2*i),*(funTable+2*i+1),*(derAnswer+2*i+1),*(t+j));
        }
        //aqui tenemos drama

        for(i=0;i<(nDimIn-j-1)*intPow(2,nDimIn-j-1);i++){
            offset = intPow(2,nDimIn-j);
            *(derAnswer+i)= linInterp(*(derAnswer+offset+2*i),*(derAnswer+offset+2*i+1),0);
        }

        j++;
    }

    //let's free some memory


    *(funTable+1) = 3; //ID method flag

    free(t); free(fl);
    return funTable;
}

double markHermite(double fpl, double fl, double fr, double fpr, double t){
double a[4]; a[0] = fl; a[1]=fpl; a[2]= -3*fl+3*fr-2*fpl-fpr; a[3]= 2*fl-2*fr+fpl+fpr;
return a[0]+t*(a[1]+t*(a[2]+t*a[3]));

}


double *get2DderTable(int *fl, double *tabConfig, double *derTable)
{
    int nDimIn = *tabConfig;
    int nDimOut = *(tabConfig +1);
    double *nBreaks = tabConfig+2;
    int nVals = getNumVals(tabConfig+2,nDimIn);

    double *answer=0; answer=(double*) calloc(sizeof(double),nDimIn*nDimOut*((int)pow(2,nDimIn))); if(answer==NULL){printf("Insuficient memory table... Closing \n");exit(1);}
    int temp=1; int j=0; int i=0; int refIndex=0;

    //this loop extracts the corner double left location of the table, so we get the refIndex.
    while(j<nDimIn-1){
        temp *=  ((int) *(nBreaks+j));j++;
        //printf("%i \n", (int) *(fl+j));
        refIndex += temp*(*(fl+j));
    }

    refIndex += (*fl); //printf("refIndex is %i \n", refIndex);
    //change the output dimension, include in refIndex
    refIndex = refIndex + (nDimOut-1)*nVals; //Read the right table
    int coordinate=0;

    //this loop fills the table with respect to refIndex. Final version.
    //It calculates the position in the table and extracts the required values
    for(i=0;i<(int) intPow(2,nDimIn);++i){
        j=0; int factor=0; temp=1; int index=0;
        while(j<nDimIn){
            //coordinate = (int) floor(i/intPow(4,j))%4; //very expensive
            coordinate =  (int) floor(i/intPow(2,j))%2; //very expensive
            factor += (coordinate * temp);
            temp *= (int) *(nBreaks+j);
            j++;
            //printf("%i ", coordinate);
        };  //printf("\n");
        index = refIndex + factor;
        //printf("%i \n", factor);
        *(answer+i) = *(derTable+index);
        //*(answer+i+4) = *(derTable+index+nVals);
    }
    return answer;
}

int *get2DCoord(double *x,double *tabConfig){
    int nDimIn = *tabConfig;
    //int nDimOut = *(tabConfig+1);
    double *nBreaks = tabConfig+2;
    int *fl = (int*) calloc(sizeof(int), (int) nDimIn); //vector of coordinates
    int i,k;
    double *min; min = (tabConfig + 2 + (int) nDimIn);
    double *max; max = (tabConfig + 2 + 2* (int) nDimIn);

    double *h = 0; h = (double*) calloc(sizeof(double),(int) nDimIn); //h is the spacing in between the data.
    for(i=0;i<nDimIn;++i){ *(h+i) = (*(max+i)- *(min+i))/( *(nBreaks+i)-1); } //get the h

    //calculate position
    for(k=0;k<nDimIn;k++){
        *(fl+k)= (int) floor((*(x+k) - *(min+k)) / (*(h+k))); //actually it is fll
       // printf("dim %i, x is %f, h %f, fl %i and t is %f \n", k+1, *(xTemp+k) , *(h+k),*(fl+k), *(t+k));
    }
    free(h);
    return fl;
}

double *getT(double *x,double *tabConfig, int *fl)
{
    int nDimIn = *tabConfig;
    //int nDimOut = *(tabConfig+1);
    double *nBreaks = tabConfig+2;
    double *t = (double*) calloc(sizeof(double), (int) nDimIn); //vector of coordinates
    int i,k;

    double *min; min = (tabConfig + 2 + nDimIn);
    double *max; max = (tabConfig + 2 + 2*nDimIn);

    double *h = 0; h = (double*) calloc(sizeof(double),nDimIn); //h is the spacing in between the data.
    for(i=0;i<nDimIn;++i){ *(h+i) = (*(max+i)- *(min+i))/(*(nBreaks+i)-1); } //get the h

    //calculate t
    for(k=0;k<nDimIn;k++){
        *(t+k)= (*(x+k)- *(min+k))/(*(h+k)) - (double) *(fl+k);
       // printf("dim %i, x is %f, h %f, fl %i and t is %f \n", k+1, *(xTemp+k) , *(h+k),*(fl+k), *(t+k));
    }
    free(h);
    return t;
}

double *get2DFun(int *fl, double *tabConfig, double *ptrTableFirst)
{   //extract hypercube from table (1D output)
    int nDimIn = *tabConfig;
    int nDimOut = *(tabConfig +1);
    double *nBreaks = tabConfig+2;
    int nVals = getNumVals(tabConfig+2,nDimIn);

    //HEISENBUG
    double *answer; answer = (double*) calloc(sizeof(double), intPow(2,nDimIn));

    if(answer==NULL){printf("Insuficient memory get2Dfun... Closing \n");exit(1);}

    int temp=1; int j=0; int i=0; int refIndex=0;

    //this loop extracts the corner double left location of the table, so we get the refIndex.
    while(j<nDimIn-1){
        temp *=  ((int) *(nBreaks+j));j++;
        //printf("%i \n", (int) *(fl+j));
        refIndex += temp*(*(fl+j));
    }

    refIndex += (*fl); //printf("refIndex is %i \n", refIndex);
    //change the output dimension, include in refIndex
    refIndex = refIndex + (nDimOut-1)*nVals; //Read the right table
    int coordinate=0;

    //this loop fills the table with respect to refIndex. Final version.
    //It calculates the position in the table and extracts the required values
    for(i=0;i<(int) intPow(2,nDimIn);++i){
        j=0; int factor=0; temp=1; int index=0;
        while(j<nDimIn){
            //coordinate = (int) floor(i/intPow(4,j))%4; //very expensive
            coordinate =  (int) floor(i/intPow(2,j))%2; //very expensive
            factor += (coordinate * temp);
            temp *= (int) *(nBreaks+j);
            j++;
            //printf("%i ", coordinate);
        };  //printf("\n");
        index = refIndex + factor;
        //printf("%i \n", factor);
        *(answer+i) = *(ptrTableFirst+index);
    }

    return answer;
}

