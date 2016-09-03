/*=========================================================================
** Title       : Interface for mexutils.c
** File        : mexutils.c
** ------------------------------------------------------------------------
** Description :
**   Interface definition header file for mexutils.c
** ------------------------------------------------------------------------ 
** Revisions   :
**   Date       Version  Author    Description
**   02-mar-06  1.0      moriborg  initial version
** -------------------------------------------------------------------------
**   (C) 2006-2010 Communication Theory Group                      
**   ETH Zurich, 8092 Zurich, Switzerland                               
**   Author: Dr. Moritz Borgmann (e-mail: moriborg@nari.ee.ethz.ch)     
** =======================================================================*/

#ifndef __MEXUTILS__
#define __MEXUTILS__

#include "mexheader.h"

typedef struct {
  long buflen;
  doublePtr r, i;
} arrBase;


void mxCopyArray(const mxArray* arr1, const mxArray* arr2);
void* getColPr(const mxArray* arr, unsigned int col);
double *getDoubleColPr(const mxArray* arr, unsigned int col);
UINT16_T *getUint16ColPr(const mxArray* arr, unsigned int col);


void getBase(mxArrayConstPtr arr, arrBase *base);
void setBase(mxArrayPtr arr, arrBase const *base);
void setBaseOffset(mxArrayPtr arr, arrBase const *base, long offset);


/* Useful accessors for array elements */
doublePtr arrGetElementPtrR(mxArrayConstPtr arr, long r, long c);
doublePtr arrGetElementPtrI(mxArrayConstPtr arr, long r, long c);

#ifdef NDEBUG
#define doubleElementR(arr, m, n) (*(mxGetPr(arr) + ((unsigned long)n)*mxGetM(arr) + ((unsigned long)m)))
#define doubleElementI(arr, m, n) (*(mxGetPi(arr) + ((unsigned long)n)*mxGetM(arr) + ((unsigned long)m)))
#else
#define doubleElementR(arr, m, n) (*arrGetElementPtrR(arr, m, n))
#define doubleElementI(arr, m, n) (*arrGetElementPtrI(arr, m, n))
#endif

#define doubleElement(arr, m, n) doubleElementR(arr, m, n)


#define uint16Element(arr, m, n) (*(((UINT16_T*)mxGetData(arr)) + ((unsigned long)n)*mxGetM(arr) + ((unsigned long)m)))

#define dumpMatrix(a) dumpMatrix_(a, #a)
void dumpMatrix_(const mxArray* arr, const char *s);

/* Debugging functions */
#define debugInt(a) mexPrintf(#a " = %d\n", a)

#define userAssert(test, msg) if(!(test)) mexErrMsgTxt(msg)
void userAssertValidArgument(const mxArray *prhs[], unsigned int ind, unsigned int m, unsigned int n, mxClassID class);

#endif
