#include <math.h>
#include "mex.h"

//  mss_sym2red.c:
//  MEX-file for [Ar,cr] = mss_sym2red(psd,A,b,c)
//
//  psd -- a 1-by-k array of positive integers indicating
//         the size of certain psd matrices (i.e. psd(i)-by-psd(i))
//
//  A   -- a m-by-P sparse double matrix.
//  b   -- a m-by-1 sparse double matrix.
//  c   -- a 1-by-P sparse double matrix
//
//  where P is the appropriate dimension to store the upper
//  triangular part of the PSD variables.
//
//  Ar -- a l-by-N sparse double matrix.
//  b  -- a l-by-1 sparse double matrix.
//  cr -- a 1-by-N sparse double matrix.
//
//  where N = prod(psd(i)^2) and l > m.


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {
  double *x,*y;              /* input x, output y      */
  mwSize m,n,q,qq,d,r,i,j,k;              /* integer counters */
  
  if(nrhs!=1) 
      mexErrMsgTxt("exactly one input required");
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) 
      mexErrMsgTxt("input must be a noncomplex double");
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  if((n!=2)||(m==0))
      mexErrMsgTxt("input must be a non-empty two-column matrix");
  x=mxGetPr(prhs[0]);
  n=m-1;
  q=0;
  i=0;
  while(i<n) {
      if((x[i+m]>=0.0)||(x[i]<=0.5))  /* not a promising set begfinning */
          i++;
      else {
          d=i+1;     /* find d: the end of neighbors <0.0 */
          while(d<m) {
              if((x[d+m]>=0.0)||(x[d]<=0.5))
                  break;
              else
                  d++;
          }
          if(x[d+m]<0.0)  /* nothing to catch in this group */
              i=d+1;
          else {
              r=d;  /* find r: the end of neighbors */
              while(r<n) {
                  if(x[r]<=0.5)
                      break;
                  else
                      r++;
              }
              q+=(d-i)*(r-d+1);
              i=r+1;
          }
      }
  }
  qq=q;
  plhs[0] = mxCreateDoubleMatrix(qq,2, mxREAL);
  y = mxGetPr(plhs[0]);
  q=0;
  i=0;
  while(i<n) {
      if((x[i+m]>=0.0)||(x[i]<=0.5))  /* not a promising set begfinning */
          i++;
      else {
          d=i+1;     /* find d: the end of neighbors <0.0 */
          while(d<m) {
              if((x[d+m]>=0.0)||(x[d]<=0.5))
                  break;
              else
                  d++;
          }
          if(x[d+m]<0.0)          /* nothing to catch in this group */
              i=d+1;
          else {
              r=d;                /* find r: the end of neighbors */
              while(r<n) {
                  if(x[r]<=0.5)
                      break;
                  else
                      r++;
              }
              for(j=i;j<d;j++)
                  for(k=d;k<=r;k++) {
                      y[q]=-x[j+m];
                      y[q+qq]=x[k+m];
                      q++;
                  }
              i=r+1;
          }
      }
  }
}
