/*=========================================================================
** Title       : Interface for matrixarith.c
** File        : matrixarith.h
** ------------------------------------------------------------------------
** Description :
**   Interface definition header file for matrixarith.c
** ------------------------------------------------------------------------ 
** Revisions   :
**   Date       Version  Author    Description
**   02-mar-06  1.0      moriborg  initial version
** -------------------------------------------------------------------------
**   (C) 2006-2010 Communication Theory Laboratory                      
**   ETH Zurich, 8092 Zurich, Switzerland                               
**   Author: Dr. Moritz Borgmann (e-mail: moriborg@nari.ee.ethz.ch)     
** =======================================================================*/

#ifndef __MATRIXARITH__
#define __MATRIXARITH__

#include "mexheader.h"

void matrixMultiply(mxArrayConstPtr R, mxArrayConstPtr A, mxArrayConstPtr B);
void matrixSubtract(mxArrayConstPtr R, mxArrayConstPtr A, mxArrayConstPtr B);
double matrixFrobeniusNormSquared(mxArrayConstPtr A);

#endif
