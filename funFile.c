#include "stdio.h" // in out
#include "stdlib.h" // for dynamic allocation
#include "math.h" // math functions
#include "time.h" // to compute time

void displayAuthor(){
    printf("Written by Anton Fadic \n");
    printf("University of Alberta \n");
    printf("Version 1.0 \n");
}

double *funTransformIn(double *inVal){
    //when this function is called, the first entry is changed as 1/T, and the others are logarithm
    *inVal=1/(*inVal);
    *(inVal+1)=log(*(inVal+1));
    *(inVal+2)=log(*(inVal+2));
    *(inVal+3)=log(*(inVal+3));
    return inVal;
}

double *invFunTransformIn(double *inVal){
    //This is invFunTransform(funTransform(x))=x
    *inVal=1/(*inVal);
    *(inVal+1)=exp(*(inVal+1));
    *(inVal+2)=exp(*(inVal+2));
    *(inVal+3)=exp(*(inVal+3));
    return inVal;
}

double *exactFun(double *inVal, int outDim, int InDim){
    //this is the exact function to compare to.
    double *temp=0; temp = malloc(sizeof(double)*InDim);
    double x,y,z;
    x=*inVal;
    y=*(inVal+1);
    // *(temp) = pow(*(inVal),2)*1; //+ (pow(*(inVal+1),2) + pow(*(inVal+2),2) ); //+ pow(*(inVal+3),2))*10000;  //+ pow(*(inVal+3),2)*1 + pow(*(inVal+4),2)*1 + pow(*(inVal+5),2)*1 ; //output dim 1
    *temp = exp((-x*x-y*y)/10);
    // *(temp+1) = 2*(*temp);
    // *(temp+2) = 3*(*temp);
    return temp;
}

void checkInBoundaries(int *fl, double *t,int nDimIn, double *nBreaks){
    //this putine changes the index *fl and *t if location is out of bounds.
    int i=0;
    for(i=0;i<nDimIn;i++){
        if(*(fl+i)<-1){*(fl+i)=-1;*(t+i)=0;}
        if(*(fl+i)>(*(nBreaks+i))-2){*(fl+i)=*(nBreaks+i)-2;*(t+i)=0;}
    }
}

int *getPosition(int *coordinate, int *nBreaks, int nDimIn){
	int i=0, *pos=0; pos = (int*) malloc(sizeof(int)*nDimIn);
	//Message("Position is: \n");
	for(i=0;i<nDimIn;i++){
		//Message("coordinate is %i nNbreaks is %i \n", *(coordinate+i), *(nBreaks+i));
		*(pos+i)=0;
		if(*(coordinate+i)<0){
			*(pos+i)=-1;
		}
		if(*(coordinate+i)+3>=*(nBreaks+i)){
			*(pos+i)=1;
		}
		//Message("%i ", *(pos+i));
	} //Message("\n");
	return pos;
}
/*
double CrInterp(double *ptrFll, double t){
    //this is the magic 4 point interpolation. This works
    double fll, fl, fr, frr;
    fll = *ptrFll; fl = *(ptrFll+1); fr = *(ptrFll+2); frr = *(ptrFll+3);
    return (1 - 3*t*t + 2*t*t*t )*fl + (3*t*t - 2*t*t*t)*fr + (t - 2*t*t + t*t*t)*(fr-fll)*0.5 + (-t*t + t*t*t)*(frr-fl)*0.5; }
*/


double CrInterp(double *ptrFll, double t, int pos){
//this is the magic 4 point interpolation
double fll, fl, fr, frr;
fll = *ptrFll; fl = *(ptrFll+1); fr = *(ptrFll+2); frr = *(ptrFll+3);
switch(pos){
	case 0 : //center
		return (1 - 3*t*t + 2*t*t*t )*fl + (3*t*t - 2*t*t*t)*fr + (t - 2*t*t + t*t*t)*(fr-fll)*0.5 + (-t*t + t*t*t)*(frr-fl)*0.5;
		break;
	case -1 : //left
		return (1 - 2*t + t*t)*fl + (2*t - t*t)*fr + (-t+t*t)*(frr-fl)*0.5;
		break;
	case 1 : //right
		return (1 - t*t)*fl + t*t*fr + (t - t*t)*(fr-fll)*0.5;
		break;
	default : //not necessary, but to prevent warning
		return (1 - 3*t*t + 2*t*t*t )*fl + (3*t*t - 2*t*t*t)*fr + (t - 2*t*t + t*t*t)*(fr-fll)*0.5 + (-t*t + t*t*t)*(frr-fl)*0.5;
	}
}


/*
double CrInterp(double *ptrFll, double t, int pos){
//this is the magic 4 point interpolation
//Use Horner's method
double fll, fl, fr, frr;
double a[4]; //polynomial coefficients.
fll = *ptrFll; fl = *(ptrFll+1); fr = *(ptrFll+2); frr = *(ptrFll+3);
//Message("%4.3f %4.3f %4.3f %4.3f t is %4.3f and pos is %i \n", fll, fl, fr, frr, t, pos);

a[0]=fl;
switch(pos){
	case 0 : //center
		a[1]=(fr-fll)*0.5; a[2]=-2.5*fl+2*fr+fll-0.5*frr; a[3]=1.5*fl-0.5*fr-1.5*fll+0.5*frr;
		//return (1 - 3*t*t + 2*t*t*t )*fl + (3*t*t - 2*t*t*t)*fr + (t - 2*t*t + t*t*t)*(fr-fll)*0.5 + (-t*t + t*t*t)*(frr-fl)*0.5;
		return a[0]+t*(a[1]+t*(a[2]+t*a[3]));
	case -1 : //left
		a[1]=-2*fl+2*fr-0.5*(frr-fl); a[2]=fl-fr+0.5*(frr-fl);
		return a[0]+t*(a[1]+t*(a[2]));
		//return (1 - 2*t + t*t)*fl + (2*t - t*t)*fr + (-t+t*t)*(frr-fl)*0.5;
	case 1 : //right
		a[1]=fr; a[2]=-fl+0.5*fr+0.5*fll;
		return a[0]+t*(a[1]+t*(a[2]));
		//return (1 - t*t)*fl + t*t*fr + (t - t*t)*(fr-fll)*0.5;
	default : //not necessary, but to prevent warning
		a[1]=(fr-fll)*0.5; a[2]=-2.5*fl+2*fr+fll-0.5*frr; a[3]=1.5*fl-0.5*fr-1.5*fll+0.5*frr;
		return a[0]+t*(a[1]+t*(a[2]+t*(a[3])));
	}
}
*/

double *readFile(int length, int nDimOut){
    FILE * fp;
    double *ptrFile = 0; ptrFile = (double*) malloc(nDimOut*length*sizeof(double));
    if(ptrFile==NULL){printf("Memory not allocated. Table could not be read. Closing"); exit(0);}
    fp = fopen ("Table.bin", "rb");
    printf("Reading table file... \n");
    fread(ptrFile,sizeof(*ptrFile),nDimOut*length,fp);
    fclose (fp);
    return ptrFile;
}

void writeTable(int length, int nDimIn, int nDimOut, double *grid){
    FILE * fp;
    int i,j;
    double *ptrTable = 0;ptrTable=(double*) malloc(nDimOut*length*sizeof(double)); if(ptrTable==NULL){printf("memory not allocated. Closing...");exit(0);}
    double *gridTemp = 0;gridTemp=(double*) malloc(sizeof(double)*nDimIn);
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
    printf("Writing table file... \n");
    fwrite(ptrTable,sizeof(*ptrTable),nDimOut*length,fp);
    free(gridTemp);
    fclose (fp);
    printf("table written successfully \n");
}

void writeTableConfig(int nDimIn, int nDimOut){
    FILE * fp;
    //int i;
    double *ptr = 0;

    ptr = (double*) malloc((3*nDimIn+2)*sizeof(double)); //nBreaks per dimension
    //write number of dimensions in
    *ptr=nDimIn;
    //write the number of dimensions out
    *(ptr+1) = nDimOut;

    //number of breaks. 0-nDimIn-1
    *(ptr+2+0) = (double) 40;
    *(ptr+2+1) = (double) 40;
    // *(ptr+2+3) = (double) 10;
    // *(ptr+2+4) = (double) 10;

    //Minimum and Maximum. 0-nDimIn-1
    *(ptr+nDimIn+2+0) = (double) 1; *(ptr+2*nDimIn+2+0) = (double) 5; //Temp
    *(ptr+nDimIn+2+1) = (double) 1; *(ptr+2*nDimIn+2+1) = (double) 5; //yNH3
    // *(ptr+nDimIn+2+2) = (double) 0.00000001; *(ptr+2*nDimIn+2+2) = (double) 0.1; //yNO
    // *(ptr+nDimIn+2+3) = (double) 0.02; *(ptr+2*nDimIn+2+3) = (double) 0.25; //yO2


    fp = fopen ("tableConfig.bin", "wb");
    fwrite(ptr,sizeof(*ptr),3*nDimIn+2,fp);
    fclose (fp);
    printf("Config file written successfully \n \n");
    free(ptr);
}

double *readTableConfig(int nDim){
    FILE * fp;
    int i;
    double *ptr2 = 0; ptr2 = (double*) malloc((3*nDim+2)*sizeof(double)); //required number of info: nDimIn, nDimOut, 3*nDim (breaks,min,max)
    if(ptr2==NULL){printf("Memory not allocated. Closing readTableConfig"); exit(0);}
    fp = fopen ("tableConfig.bin", "rb");
    printf("Reading config file... \n");
    fread(ptr2,1,(3*nDim+2)*sizeof(*ptr2),fp);
    fclose (fp);
    printf("Config file read successfully \n \n");
    printf("Number of input dimensions is %i \n", (int) *ptr2);
    printf("Number of output dimensions is %i \n", (int) *(ptr2+1));
    for(i=0;i<nDim;++i){
        printf("Dimension %i, breaks = %i, Min = %f , Max = %f \n", i+1,(int) *(ptr2+i+2), *(ptr2+nDim+i+2), *(ptr2+2*nDim+i+2));
    }
    return ptr2;
}

double *writeGrid(double *ptrFirstVal, int nDimIn, int nVals){
    //This function is in its final version. Please do not touch. This function fails for more than 6 dimensions
    //because of segmentation fault. Working on this. Possible workaraounds, fragmentate de memory.
    int i,j;
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
    return grid;
}

int getNumVals(double *ptrFirstVal, int nDim){
    // *ptrFirstVal is pointing to the number of breaks
    int i;
    int nVals=1;
    double dimsArray[nDim];
    for(i=0;i<nDim;++i){
       dimsArray[i] = *(ptrFirstVal+i);
       nVals = nVals*dimsArray[i];
    }
    printf("number of grid values per output dimension is: %i \n", nVals);
    return nVals;
}

void saveGrid(double *grid, int nDimIn, int nVals){
    FILE *fp;
    fp = fopen("grid.bin","wb");
    fwrite(grid,sizeof(double),nVals*nDimIn,fp);
    fclose(fp);
}

double *interpolate(double *x, int nDimIn, int nDimOut, double *ptrFirstVal, double *ptrTableFirst, int nVals){
    //this function interpolates values one by one.
    int i,j,k;
    double *nBreaks; nBreaks = ptrFirstVal+2;
    double *min; min = (ptrFirstVal + 2 + nDimIn);
    double *max; max = (ptrFirstVal + 2 + 2*nDimIn);

    //apply transformations
    /*
    min = funTransformIn(min);
    min = invFunTransformIn(min);
    max = funTransformIn(max);
    max = invFunTransformIn(max);
    x = funTransformIn(x);
    x = invFunTransformIn(x);
    */
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
    pos = getPosition(fl, (int*) nBreaks, nDimIn); //get the position if it lies of the borders for interpolation order
    checkInBoundaries(fl,t,nDimIn,nBreaks);       //adjust position if it lies outside borders. Safety measure to avoid extrapolation.

    //extract hypercube from table (1D output)
    double *answer=0; answer=(double*) malloc(sizeof(double)*nDimOut*((int)pow(4,nDimIn))); if(answer==NULL){printf("Insuficient memory interpolate... Closing \n");exit(1);}
    int refIndex=0; int temp=1; j=0;

    //this loop extracts the corner double left location of the table, so we get the refIndex.
    while(j<nDimIn-1){
        temp *=  ((int) *(nBreaks+j));j++;
        //printf("%i \n", (int) *(fl+j));
        refIndex += temp*(*(fl+j));
    }

    refIndex += (*fl); //printf("refIndex is %i \n", refIndex);
    //change the output dimension, include in refIndex
    refIndex = refIndex + (nDimOut-1)*nVals; //Read the right table
    int index=0; int coordinate=0; int factor;

    //this loop fills the table with respect to refIndex. Final version.
    //this is the most expensive function of the code!!!
    //It searches the table and extracts the required values
    //Takes ~90% of time.
    for(i=0;i<(int) pow(4,nDimIn);++i){
        j=0; factor=0; temp=1;
        while(j<nDimIn){
            coordinate = (int) floor(i/pow(4,j))%4; //do not touch
            factor += (coordinate * temp);
            temp *= ((int) *(nBreaks+j));
            j++;
            //printf("%i ", coordinate);
        };  //printf("\n");
        index = refIndex + factor;
        //printf("%i \n", factor);
        *(answer+i) = *(ptrTableFirst+index);
    }

    //this loop creates the contraction to interpolate and the value is
    //kept in the first position of the memory. This only takes ~10% of time.
    j=0;

    while(j<nDimIn){
        //printf("t is %f \n", *(t+j));
        i=0;
        while(i<(int) pow(4,nDimIn-j-1)){
             *(answer+i) = CrInterp(answer+4*i,*(t+j),*(pos+j));
            // *(answer+i) = CrInterp(answer+4*i,*(t+j));
            //printf("int Calc %f \n", *(answer+i));
            ++i;
        }
        j++;
    }
    return answer;
    // *answer = exp((*answer)/1000000);
    //undo transformation;
    //x = invFunTransformIn(x);
    //min = invFunTransformIn(min);
    //max = invFunTransformIn(max);


    //printf("answer is %f \n",*answer);
    //printf("exact  is: %f \n", *(exactFun(x,nDimOut,nDimIn)+nDimOut-1));
    //printf("error is %f \% \n", 100-100*(*answer/(*(exactFun(x,nDimOut,nDimIn)+nDimOut-1))) );
   // free(answer); free(h); free(t); free(xTemp); free(fl);
}
