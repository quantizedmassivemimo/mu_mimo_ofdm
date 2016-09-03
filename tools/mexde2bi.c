/*=========================================================================
** Title       : Decimal to Binary Conversion
** File        : mexde2bi.c
** ------------------------------------------------------------------------
** Description :
**   Replacement for de2bi, expects a column vector, and the number
**   of bits. Returns a matrix, with the number of rows equal to the
**   input vector and the bits in the columns
** ------------------------------------------------------------------------ 
** Revisions   :
**   Date       Version  Author  Description
**   23-nov-05  1.0      mimo    file created
** -------------------------------------------------------------------------
**   (C) 2006-2010 Communication Theory Laboratory                      
**   ETH Zurich, 8092 Zurich, Switzerland                               
**   Author: Dr. Christoph Studer (e-mail: studerc@nari.ee.ethz.ch)     
** =======================================================================*/

#ifndef __MEXHEADER__
#define __MEXHEADER__

#include <mex.h>
#include <math.h>
#include <string.h>

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define NIL 0L
#define nil 0L

typedef mxArray * mxArrayPtr;
typedef const mxArray * mxArrayConstPtr;
typedef double * doublePtr;

#endif

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A;
    unsigned int rowsA,colsA, maxdimA;
    unsigned int a;
    int NrOfBits,bit;
    double *R;
    int i;
    
    /* Get the number of rows in A and the data pointers */
    A = mxGetPr(prhs[0]);
    rowsA =  mxGetM(prhs[0]);
    colsA =  mxGetN(prhs[0]);
    maxdimA = (rowsA > colsA) ? rowsA : colsA;
    NrOfBits=*mxGetPr(prhs[1]);
    /* Create the result array */
    plhs[0] = mxCreateDoubleMatrix(maxdimA, NrOfBits, mxREAL);
    R = mxGetPr(plhs[0]);

    
    for(i=0;i<maxdimA;i++) {
        /* de2bi */
        a=(int)A[i];
        for(bit=0;bit<NrOfBits;bit++) {
            if(a & 0x1) {
                *(R + i + maxdimA*bit) = 1.0;
            } else {
                *(R + i + maxdimA*bit) = 0.0;
            }
            a = a >> 1;
        }
    }
    
}
