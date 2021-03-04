/* edcmp: compute elastic Green's functions for a layered half space.
 *   Sylvain Barbot
 *   sbarbot@ntu.edu.sg
 *   Earth Observatory of Singapore
*/

//#include <ctype.h>
#include "mex.h"
#include "matrix.h"
//using namespace std;

extern"C" {
void edcmp_(double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, double *, double *, double *, double *, double *, int *);
}

void mexFunction(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[]) {

  // function [ux, uy, uz] = layered_wang(slip,xs,ys,zs,L,W,strike,dip,rake,xr,yr)
  double * x, *y;
  double * u1, *u2, *u3;
  size_t nxx, nxy, ny;
  int nrec,err;

  if (9>nrhs)
     mexErrMsgIdAndTxt("MATLAB:edcmp:minrhs",
                       "not enough input arguments.");

  // test input
  double slip = mxGetScalar(prhs[0]);
  double x1 = mxGetScalar(prhs[1]);
  double x2 = mxGetScalar(prhs[2]);
  double x3 = mxGetScalar(prhs[3]);
  double length = mxGetScalar(prhs[4]);
  double width = mxGetScalar(prhs[5]);
  double strike = mxGetScalar(prhs[6]);
  double dip = mxGetScalar(prhs[7]);
  double rake = mxGetScalar(prhs[8]);

  mexPrintf("slip=%f\nx1=%f\nx2=%f\nx3=%f\nL=%f\nW=%f\nstrike=%f\ndip=%f\nrake=%f\n",
            slip,x1,x2,x3,length,width,strike,dip,rake);  

  if (3>nlhs) 
     mexErrMsgIdAndTxt("MATLAB:edcmp:minlhs",
                       "not enough output arguments.");

  /* get the length of each input vector */
  nxx = mxGetM(prhs[9]);
  ny = mxGetN(prhs[9]);

  if (1 != ny)
     mexErrMsgIdAndTxt("MATLAB:edcmp:dimension",
                       "x coordinates should be vectors.");

  nxy = mxGetM(prhs[10]);
  ny = mxGetN(prhs[10]);

  if (1 != ny)
     mexErrMsgIdAndTxt("MATLAB:edcmp:dimension",
                       "y coordinates should be vectors.");

  if (nxx != nxy)
     mexErrMsgIdAndTxt("MATLAB:edcmp:dimension",
                       "x and y coordinates should have the same dimension.");

  plhs[0]=mxCreateDoubleMatrix((mwSize)nxx,(mwSize)ny,mxREAL);
  plhs[1]=mxCreateDoubleMatrix((mwSize)nxx,(mwSize)ny,mxREAL);
  plhs[2]=mxCreateDoubleMatrix((mwSize)nxx,(mwSize)ny,mxREAL);

  u1=mxGetPr(plhs[0]);
  u2=mxGetPr(plhs[1]);
  u3=mxGetPr(plhs[2]);

  x=mxGetPr(plhs[9]);
  y=mxGetPr(plhs[10]);

  nrec=(int)nxx;
  edcmp_(&slip,&x1,&x2,&x3,&length,&width,&strike,&dip,&rake,&nrec,x,y,u1,u2,u3,&err);
}

