/*=========================================================================
** Title       : Useful Definitions for Matlab MEX
** File        : mexheader.c
** ------------------------------------------------------------------------
** Description :
**   C header file with several useful definitions for MATLAB MEX 
**   programming.
** ------------------------------------------------------------------------ 
** Revisions   :
**   Date       Version  Author    Description
**   02-mar-06  1.0      moriborg  initial version
** -------------------------------------------------------------------------
**   (C) 2006-2010 Communication Theory Group                      
**   ETH Zurich, 8092 Zurich, Switzerland                               
**   Author: Dr. Moritz Borgmann (e-mail: moriborg@nari.ee.ethz.ch)     
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
