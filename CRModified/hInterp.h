double CrInterp(double *ptrFll, double t); //interpolates

double *readFile(int length, int nDimOut); //read the table

void writeTable(int length, int nDimOut, int nDimIn, double *grid); // writes the table

void writeTableConfig(int nDim, int nDimOut, int brX, int brY); //write config file of table.

double *readTableConfig(int nDim); //read the table config

double *writeGrid(double *ptrFirstVal); // writes the grid

int getNumVals(double *ptrFirstVal, int nDim); // get the number of values from the number of dimentions

void saveGrid(double *grid, int nDimIn, int nVals); // stores the grid in Disk

double *interpolate(double *x, double *ptrFirstVal, double *ptrTableFirst); //the most important function in the world

double *funTransformIn(double *inVal);

double *invFunTransformIn(double *inVal);

void checkInBoundaries(int *fl, double *t,int nDimIn);

void displayMethod(double *answer);

int *getPosition(int *coordinate, double *nBreaks, int nDimIn);

double *calcDeriv(int nDimIn, int nDimOut, double *answer);

double HermiteInterp(double fll, double fl, double fr, double frr, double t, int pos);

double *getFunVals(int nDimIn, double* answer);

int intPow(int a, int b);

void displayAuthor();

double *interpHermite(double *x, double *tabConfig, double *ptrTableFirst,  double *derTable);

double linInterp(double fl, double fr, double t);

double *get2DderTable(int *fl, double *tabConfig, double *ptrTableFirst);

int *get2DCoord(double *x,double *tabConfig);

double *getT(double *x,double *tabConfig, int *fl);

double *get2DFun(int *fl, double *tabConfig, double *ptrTableFirst);

double *getAnswer(int nDimIn, double *ptrTableFirst, int *fl, int nDimOut, double *nBreaks, int nVals);

double *getDerTable(double *tabConfig, double *grid);

double *PolyMultiDim(int nDimIn, int nDimOut, double *answer, double *t, int *pos);

double markHermite(double fpl, double fl, double fr, double fpr, double t);
