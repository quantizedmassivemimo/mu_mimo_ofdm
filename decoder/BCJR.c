/*=========================================================================
** Title       : Max-log BCJR Implementation
** File        : BCJR.c
** ------------------------------------------------------------------------
** Description :
**   Fast implementation of a max-log soft-in soft-out MAP decoder 
**   [Bahl et.al.] aka. forward-backward algorithm. Decoder has two 
**   outputs: a posteriori LLRs for coded bits and sliced output estimates 
**   (bits). 
** ------------------------------------------------------------------------ 
** Revisions   :
**   Date       Version  Author  Description
**   17-aug-06  1.3      studer  some bugfixes/cleanup
**   12-jul-06  1.2      studer  minor speed-ups
**   01-jun-06  1.1      studer  a posteriori output and minor speedups
**   01-may-06  1.0      studer  initial version adapted from softviterbi
** -------------------------------------------------------------------------
**   (C) 2006-2010 Communication Theory Group                      
**   ETH Zurich, 8092 Zurich, Switzerland                               
**   Author: Dr. Christoph Studer (e-mail: studerc@nari.ee.ethz.ch)     
** =======================================================================*/

#include "mexheader.h"
#include "mexutils.h"
#include "matrixarith.h"
#include "math.h"

unsigned int log2floor(unsigned int x) 
{
    unsigned int result;
    result=0;    
    while((1 << result) <= x) result++;
    result--;
    return result;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int 	          k, K, i, j, q;
  unsigned int    p;
  
  mxArrayConstPtr mtrellisArr, tempArr;
  mxArrayPtr      outSymbolsArr, nextStatesArr;
  mxArrayConstPtr LLRVector;
  mxArrayPtr      BMVector;
  
  unsigned int    numOutputSymbols, numInputSymbols;
  unsigned int    numOutputBits, numInputBits;
  unsigned int    numStates;
  unsigned int    numBits, numSymbols;

  mxArrayPtr	  AlphaMetric, BetaMetric, oldBetaMetric;
  mxArrayPtr      xhatk, xhatk0, xhatk1;
  mxArrayPtr      llr_0, llr_1, llr_d;

  /* Examine input and output arguments. */
  userAssert(nrhs == 2, "BCJR.c requires 2 input arguments!");    
  userAssert(nlhs == 2, "BCJR.c requires 2 output arguments!");
  
  /* Read out mtrellis structure into variables */
  mtrellisArr = prhs[0];
  
  tempArr = mxGetField(mtrellisArr, 0, "numInputSymbols");
  mxAssert(tempArr != NIL, "Missing field 'numInputSymbols' in mtrellis parameter");
  numInputSymbols = (unsigned int)(*mxGetPr(tempArr));
  numInputBits = (int)log2floor(numInputSymbols);

  tempArr = mxGetField(mtrellisArr, 0, "numOutputSymbols");
  mxAssert(tempArr != NIL, "Missing field 'numOutputSymbols' in mtrellis parameter");
  numOutputSymbols = (unsigned int)*mxGetPr(tempArr);
  numOutputBits = (int)log2floor(numOutputSymbols);

  tempArr = mxGetField(mtrellisArr, 0, "numStates");
  mxAssert(tempArr != NIL, "Missing field 'numStates' in mtrellis parameter");
  numStates = (unsigned int)*mxGetPr(tempArr);

  nextStatesArr = mxGetField(mtrellisArr, 0, "nextStates");
  userAssert(nextStatesArr != NIL, "Missing field 'nextStates' in mtrellis parameter");
  
  outSymbolsArr = mxGetField(mtrellisArr, 0, "outputs");
  userAssert(outSymbolsArr != NIL, "Missing field 'outputs' in mtrellis parameter");
  
  /* Get received bits */
  LLRVector = prhs[1];
  if(mxGetM(LLRVector) > mxGetN(LLRVector)) numBits = mxGetM(LLRVector); 
    else numBits = mxGetN(LLRVector); 
  numSymbols = numBits/numOutputBits;

  /* -------------------------------- */
  
  /* detected systematic message */
  xhatk0 = mxCreateDoubleMatrix(1, numInputBits*numSymbols, mxREAL);
  xhatk1 = mxCreateDoubleMatrix(1, numInputBits*numSymbols, mxREAL);
  xhatk  = mxCreateDoubleMatrix(1, numInputBits*numSymbols, mxREAL);
  
  /* a posteriori llr output */
  llr_0 = mxCreateDoubleMatrix(1,numBits,mxREAL);
  llr_1 = mxCreateDoubleMatrix(1,numBits,mxREAL);
  llr_d = mxCreateDoubleMatrix(1,numBits,mxREAL);
     
  /* Branch Metric Vector (precomputed) */
  BMVector = mxCreateDoubleMatrix(numOutputSymbols, 1, mxREAL);

  /* ------------------------------------------------------------------- */
  /* Alpha-Metric (Forward Iteration)                                    */
  /* ------------------------------------------------------------------- */
  AlphaMetric    = mxCreateDoubleMatrix(numStates, (int)(numSymbols+1), mxREAL); /* initializes with zero */
  /* By definition, we start in state 0. In order to make this clear
  ** to the BCJR, we severely discriminate against all other states... */
  for(i = 0; i<numSymbols; i+=1)
      for(j = 0;j<numStates; j++)
          doubleElement(AlphaMetric, i, j) = -DBL_MAX;
  
  doubleElement(AlphaMetric, 0, 0) = 100000.0;
  
  for (i = 0; i<numSymbols; i+=1) {
    int tmp,q,b;
    /* precompute the Branch Metric for each output Symbol, based on the LLRs */
    for(q=0;q<numOutputSymbols;q++) {
      doubleElementR(BMVector,q,0) = (double) 0.0;
      tmp=q;
      for(b=numOutputBits-1;b>=0;b--) {
        if(tmp>=(1 << b)) {
          tmp = tmp - (1 << b);
          doubleElementR(BMVector,q,0) += doubleElementR(LLRVector,(i*numOutputBits+numOutputBits-1-b),0);   /* first bit is MSB */
        } else {
          doubleElementR(BMVector,q,0) -= doubleElementR(LLRVector,(i*numOutputBits+numOutputBits-1-b),0);   /* first bit is MSB */
        }
      }
      /* adjust to max-log metric */
      doubleElementR(BMVector,q,0) = 0.5*doubleElementR(BMVector,q,0);
    }
    
    /* ------------------ carry out one trellis step ------------------------- */
    for (j = 0;j<numStates; j++) {                            /* for each to-state */
      double prevAlpha, curMetric;
      prevAlpha = doubleElementR(AlphaMetric, j, i);      
      /* for each branch */
      for (p = 0; p <numInputSymbols; p++) {                  /* for each branch */
	    int nextState;	    
        curMetric = prevAlpha + doubleElementR(BMVector,doubleElementR(outSymbolsArr,j,p),0);
        nextState = doubleElementR(nextStatesArr,j,p);
        if(curMetric > doubleElementR(AlphaMetric,nextState,i+1)) doubleElementR(AlphaMetric, nextState, i+1) = curMetric;
      } /* for p */
    } /* for j */
  } /* for i */
  
  /* ------------------------------------------------------------------- */
  /* Beta-Metric (Backward Iteration)                                    */
  /* ------------------------------------------------------------------- */
  BetaMetric    = mxCreateDoubleMatrix(numStates, 1, mxREAL);
  oldBetaMetric = mxCreateDoubleMatrix(numStates, 1, mxREAL);
  for(j = 0;j<numStates; j++) doubleElement(oldBetaMetric, j, 0) = -DBL_MAX;
  doubleElement(oldBetaMetric, 0, 0) = 100000.0; /* Termination */
  for (i = numSymbols; i>0; i-=1) {
    int tmp,q,b;
    /* precompute the Branch Metric for each output Symbol, based on the LLRs */
    for(q=0;q<numOutputSymbols;q++) {
      doubleElementR(BMVector,q,0) = (double)0.0;
      tmp=q;
      for(b=numOutputBits-1;b>=0;b--) {
        if(tmp>=(1 << b)) {
          tmp = tmp - (1 << b);
          doubleElementR(BMVector,q,0) += doubleElementR(LLRVector,(i*numOutputBits-1-b),0);   /* first bit is MSB */
        } else {
          doubleElementR(BMVector,q,0) -= doubleElementR(LLRVector,(i*numOutputBits-1-b),0);   /* first bit is MSB */
        }
      }
      /* adjust to max-log metric */
      doubleElementR(BMVector,q,0) = 0.5*doubleElementR(BMVector,q,0);
    }
      
    /* ------------------ carry out one trellis step ------------------------- */
    for (j = 0;j<numStates; j++) {                            /* for each state */
      double metricInc, curMetric;
      double outputMetric; 
      int bit; 
      /* for each branch */
      for (p = 0; p <numInputSymbols; p++)  { /* for each branch */
        /* ---- Beta ------------- */
        int nextState;
        metricInc = doubleElementR(BMVector,doubleElementR(outSymbolsArr,j,p),0);	
        nextState = doubleElementR(nextStatesArr,j,p);
        curMetric = doubleElementR(oldBetaMetric,(int)nextState,0) + metricInc;
        if ((curMetric > doubleElementR(BetaMetric,j,0)) || (p==0)) 
          doubleElementR(BetaMetric, j, 0) = curMetric;

        /* ---- output metric of current state ---- */
        outputMetric = curMetric + doubleElementR(AlphaMetric, j, i-1); 
          
        /* ---- LLR output of input bits ---- */
        for (b=numInputBits-1;b>=0;b--) 
        {
          int BitIdx;
          BitIdx = numInputBits*i-1-b;
          if (p & (1<<b)) 
          { /* bit is +1 */
            if (outputMetric>doubleElementR(xhatk1,BitIdx,0)) doubleElementR(xhatk1,BitIdx,0) = outputMetric;
          } else 
          { /* bit is -1 */          
            if (outputMetric>doubleElementR(xhatk0,BitIdx,0)) doubleElementR(xhatk0,BitIdx,0) = outputMetric;
          }
        }
                                  
        /* ---- a posteriori LLR output calculation ---- */                
        bit = (int) doubleElementR(outSymbolsArr,j,p);        
        for (b=numOutputBits-1;b>=0;b--)
        {
          int BitIdx;
          BitIdx = i*numOutputBits-1-b;          
          if (bit>=(1<<b))
          { /* bit is +1 */
            bit = bit - (1<<b);
            if (outputMetric>doubleElementR(llr_1,BitIdx,0)) doubleElementR(llr_1,BitIdx,0) = outputMetric;
          } else
          { /* bit is -1 */
            if (outputMetric>doubleElementR(llr_0,BitIdx,0)) doubleElementR(llr_0,BitIdx,0) = outputMetric;
          }
        }                  
      
      } /* for p input symbols */
    } /* for j states */
    
    /* copy the beta state metric */
    for(q=0;q<numStates;q++) 
      doubleElementR(oldBetaMetric,q,0)=doubleElementR(BetaMetric,q,0);
    
  } /* for i branches */

  /* Compute the sliced systematic output bits */
  for (i=0;i<numInputBits*numSymbols;i++)  
    if (doubleElementR(xhatk1,i,0)>doubleElementR(xhatk0,i,0))
      doubleElementR(xhatk,i,0) = 1.0;
    else
      doubleElementR(xhatk,i,0) = 0.0;       
          
  /* Compute a posteriori information [LLRs] */
  for (i=0;i<numBits;i++) 
    doubleElementR(llr_d,i,0) = doubleElementR(llr_1,i,0)-doubleElementR(llr_0,i,0);
  
  /* output */   
  plhs[1] = xhatk; /* -- systematic bits output -- */
  plhs[0] = llr_d; /* -- a posteriori output -- */ 
  
  /* free memory */  
  mxDestroyArray(BetaMetric);
  mxDestroyArray(oldBetaMetric);
  mxDestroyArray(xhatk0);
  mxDestroyArray(xhatk1);  
  mxDestroyArray(llr_0);
  mxDestroyArray(llr_1);    
  
}
