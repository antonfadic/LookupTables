#include "stdio.h" // in out
#include "math.h" // math functions
#include "exactFun.h" // include exact function
#include "stdlib.h" // for dynamic allocation

void checkInBoundaries(int *fl, double *t,int nDimIn, double *nBreaks){
    //this routine changes the index *fl and *t if location is out of bounds.
    int i=0;
    for(i=0;i<nDimIn;i++){
        if(*(fl+i)<-1){*(fl+i)=-1;*(t+i)=0;}
        if(*(fl+i)>(*(nBreaks+i))-2){*(fl+i)=*(nBreaks+i)-2;*(t+i)=0;}
    }
}

int *getPosition(int *coordinate, double *nBreaks, int nDimIn){
	int i=0, *pos=0; pos = (int*) calloc(sizeof(int),nDimIn);
	//Message("Position is: \n");
	for(i=0;i<nDimIn;i++){
		//printf("coordinate is %i nNbreaks is %i \n", *(coordinate+i), *(nBreaks+i));
		*(pos+i)=0;
		if(*(coordinate+i)+1<=0){
			*(pos+i)=-1;
		}
		if(*(coordinate+i)+3>=(int) *(nBreaks+i)){
			*(pos+i)=1;
		}
		//Message("%i ", *(pos+i));
	} //Message("\n");
	return pos;
}

double HermiteInterp(double fpl, double fl, double fr, double fpr, double t, int pos){

    double a[4]; a[0] = fl; a[1]=fpl; a[2]= -3*fl+3*fr-2*fpl-fpr; a[3]= 2*fl-2*fr+fpl+fpr;
    pos=0;
    switch(pos){
        case 0 : //center
		//return a[0] + a[1]*t + a[2] *t*t + a[3] *t*t*t;
		return a[0]+t*(a[1]+t*(a[2]+t*a[3]));
        break;
	case -1 : //left
	    a[1]= 2*fr - 2*fl - fpr; a[2]=fl - fr + fpr;
		return a[0]+t*(a[1]+t*a[2]);
		break;
	case 1 : //right
	    a[2]= fr - fl - fpl;
		return a[0]+t*(a[1]+t*a[2]);
		break;
	default : //not necessary, but to prevent warning
		return a[0] + a[1]*t + a[2] *t*t + a[3] *t*t*t;
	}
}

double CrInterp(double *ptrFll, double t, int pos){
//this is the magic 4 point interpolation
double fll, fl, fr, frr, fpl, fpr;
fll = *ptrFll; fl = *(ptrFll+1); fr = *(ptrFll+2); frr = *(ptrFll+3);
fpl = 0.5*(fr-fll); fpr=0.5*(frr-fl);
double a[4]; a[0] = fl; a[1]=fpl; a[2]= -3*fl+3*fr-2*fpl-fpr; a[3]= 2*fl-2*fr+fpl+fpr;

switch(pos){
	case 0 : //center
		//return (1 - 3*t*t + 2*t*t*t )*fl + (3*t*t - 2*t*t*t)*fr + (t - 2*t*t + t*t*t)*(fr-fll)*0.5 + (-t*t + t*t*t)*(frr-fl)*0.5;
		return a[0]+t*(a[1]+t*(a[2]+t*a[3]));
		break;
	case -1 : //left
        a[1]= 2*fr - 2*fl - fpr; a[2]=fl - fr + fpr;
		return a[0]+t*(a[1]+t*a[2]);
		break;
	case 1 : //right
		a[2]= fr - fl - fpl;
		return a[0]+t*(a[1]+t*a[2]);
		break;
	default : //not necessary, but to prevent warning
		return (1 - 3*t*t + 2*t*t*t )*fl + (3*t*t - 2*t*t*t)*fr + (t - 2*t*t + t*t*t)*(fr-fll)*0.5 + (-t*t + t*t*t)*(frr-fl)*0.5;
	}
}

double CrInterpAlt(double *ptrFll, double t, int pos){
//this is the magic 4 point interpolation. This is even more magic.
double fll, fl, fr, frr, fpl, fpr;
fll = *ptrFll; fl = *(ptrFll+1); fr = *(ptrFll+2); frr = *(ptrFll+3);
fpl = -fll/3-fl/2+fr-frr/6; fpr=fll/6-fl+fr/2+frr/3;
double a[4]; a[0] = fl; a[1]=fpl; a[2]= -3*fl+3*fr-2*fpl-fpr; a[3]= 2*fl-2*fr+fpl+fpr;

switch(pos){
	case 0 : //center
		//return (1 - 3*t*t + 2*t*t*t )*fl + (3*t*t - 2*t*t*t)*fr + (t - 2*t*t + t*t*t)*(fr-fll)*0.5 + (-t*t + t*t*t)*(frr-fl)*0.5;
		return a[0]+t*(a[1]+t*(a[2]+t*a[3]));
		break;
	case -1 : //left
	    fpr=(frr-fl)*0.5;
        a[1]= 2*fr - 2*fl - fpr; a[2]=fl - fr + fpr;
		return a[0]+t*(a[1]+t*a[2]);
		break;
	case 1 : //right
		a[2]= fr - fl - fpl;
		fpl = (fr-fll)/2;
		return a[0]+t*(a[1]+t*a[2]);
		break;
	default : //not necessary, but to prevent warning
		return (1 - 3*t*t + 2*t*t*t )*fl + (3*t*t - 2*t*t*t)*fr + (t - 2*t*t + t*t*t)*(fr-fll)*0.5 + (-t*t + t*t*t)*(frr-fl)*0.5;
	}
}

int intPow(int a, int b){
    int i, result;
    result=1;
    for(i=0;i<b;i++){
        result *= a;
    }
    return result;
}

double *calcDeriv(int nDimIn, int nDimOut, double *answer) {
    int i,j,k, nVals, coordinate, index, factor, temp, dummy, coordinateL, factorL;
    k=0; nVals=1;
    //for(i=0;i<nDimIn;i++){nVals*=2;}; //number of values per dimension
    nVals = intPow(2,nDimIn);
    double *derGrid =0; derGrid = (double*) calloc(sizeof(double),nDimIn*nVals);
    while(k<nDimIn){
        for(i=0;i<nVals;++i){
            j=0; factor=0; temp=1; dummy=0; factorL=0;
            while(j<nDimIn){
                if(k==j){dummy=1;}; if(k!=j){dummy=0;}
                coordinate = (int) floor(i/intPow(2,j))%2 + 1+ 1*dummy; //do not touch
                coordinateL = (int) floor(i/intPow(2,j))%2 + 1- 1*dummy; //do not touch
                factor += (coordinate * temp);
                factorL += (coordinateL * temp);
                temp *= 4;
                j++;
                //printf("%i ", coordinateL);
            };  // printf("\n");
            index = factor;
            //printf("%i %i \n", factor, factorL);
            *(derGrid+i+k*nVals) = *(answer+index)*0.5 - *(answer+factorL)*0.5;
        }k++;
    }
    return derGrid;
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
    //printf("number of grid values per output dimension is: %i \n", nVals);
    return nVals;
}

double *getFunVals(int nDimIn, double* answer){
//this calculates only the 2^D hypercube
double *funTable =0; funTable = (double*) calloc(sizeof(double),intPow(2,nDimIn));
int i,j, coordinate, index, factor, temp;
   for(i=0;i<(int) intPow(2,nDimIn);++i){
        j=0; factor=0; temp=1;
        while(j<nDimIn){
            coordinate = (int) floor(i/intPow(2,j))%2 + 1; //do not touch
            factor += (coordinate * temp);
            temp *= 4;
            j++;
            //printf("%i ", coordinate);
        };  //printf("\n");
        index = factor;
        //printf("%i \n", factor);
        *(funTable+i) = *(answer+index);
    }
    return funTable;
}

double linInterp(double fl, double fr, double t){
    return fl + t*(fr-fl);
}

int intFloor(double x){
    int i= (int) x;
    return i - ( i > x );
}

double *CRMultiDim(int nDimIn, int nDimOut, double *answer, double *t, int *pos){
    int i,j=0;
    while(j<nDimIn){
        //printf("t is %f \n", *(t+j));
        i=0;
        while(i<(int) intPow(4,nDimIn-j-1)){
             *(answer+i) = CrInterp(answer+4*i,*(t+j),*(pos+j));
            //printf("int Calc %f \n", *(answer+i));
            ++i;
        }
        j++;
    }
    //let's free some memory
    answer = (double*) realloc(answer,sizeof(double)*nDimOut*2);
    *(answer+1) = 1.0; //flag for chosen method
    return answer;
}

double *PolyMultiDim(int nDimIn, int nDimOut, double *answer, double *t, int *pos){
    int i,j=0;
    while(j<nDimIn){
        //printf("t is %f \n", *(t+j));
        i=0;
        while(i<(int) intPow(4,nDimIn-j-1)){
             *(answer+i) = CrInterpAlt(answer+4*i,*(t+j),*(pos+j));
            //printf("int Calc %f \n", *(answer+i));
            ++i;
        }
        j++;
    }
    //let's free some memory
    answer = (double*) realloc(answer,sizeof(double)*nDimOut*2);
    *(answer+1) = 2.0;
    return answer;
}

double *CombiMultiDim(int nDimIn, int nDimOut, double *answer, double *t, int *pos){
    //combined method
    double *funTable=0; funTable = getFunVals(nDimIn,answer); //this gets the 2D function of values
    double *derGrid=0; derGrid = calcDeriv(nDimIn,nDimOut,answer); //calculate the derivatives
    int i,j,offset; i=0;j=0;offset=0;

    while(j<nDimIn){
        for(i=0;i<intPow(2,nDimIn-j-1);i++){
            *(funTable+i) = HermiteInterp(*(derGrid+2*i),*(funTable+2*i),*(funTable+2*i+1),*(derGrid+2*i+1),*(t+j),*(pos+j));
        }
        for(i=0;i<(nDimIn-j-1)*intPow(2,nDimIn-j-1);i++){
            offset = intPow(2,nDimIn-j);
            *(derGrid+i)= linInterp(*(derGrid+offset+2*i),*(derGrid+offset+2*i+1),*(t+j));
        }
        j++;
    }
    //let's free some memory
    free(answer); free(derGrid);
    *(funTable+1) = 0;
    return funTable;
}

double *getAnswer(int nDimIn, double *ptrTableFirst, int *fl, int nDimOut, double *nBreaks, int nVals){
   //extract hypercube from table (1D output)
    double *answer=0; answer=(double*) malloc(sizeof(double)*nDimOut*((int)pow(4,nDimIn))); if(answer==NULL){printf("Insuficient memory interpolate... Closing \n");exit(1);}
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
    int index=0; int coordinate=0; int factor;

    //this loop fills the table with respect to refIndex. Final version.
    //It calculates the position in the table and extracts the required values
    for(i=0;i<(int) intPow(4,nDimIn);++i){
        j=0; factor=0; temp=1;
        while(j<nDimIn){
            //coordinate = (int) floor(i/intPow(4,j))%4; //very expensive
            coordinate = intFloor(i/intPow(4,j))%4; //very expensive
            factor += (coordinate * temp);
            temp *= ((int) *(nBreaks+j));
            j++;
            //printf("%i ", coordinate);
        };  //printf("\n");
        index = refIndex + factor;
        //printf("%i \n", factor);
        *(answer+i) = *(ptrTableFirst+index);
    }
    return answer;
}

void displayAuthor(){
    printf("Written by Anton Fadic \n");
    printf("University of Alberta \n");
    printf("Modified version, combined Hermite-CR\n");
    printf("Version 1.0 beta \n");
}

void displayMethod(double *answer){
    int fCase;
    fCase = (int) *answer;
    printf("\n");
    switch(fCase) {
    case 0:
        printf("Combined Multi dimensional interpolation \n");
        break;
    case 1:
        printf("Pure CR Interpolation \n");
        break;
    case 2:
        printf("4 Point poly spline \n");
        break;
    case 3:
        printf("Markus method \n");
        break;
    }
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

