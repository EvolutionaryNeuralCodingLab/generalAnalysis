/*[C]=crossCorrRaster(t,ic,maxLag)
 * t must be rounded.
 * cross correlation of all spike trains.
 * To compile: mex -largeArrayDims crossCorrRaster.c
 * Format:
 * [C]=crossCorrRaster(t,ic,lag);
 * t is the time vector
 * ic is the index channel
 * lag is the desired lag
 * and C is the calculated cross correlation matrix
 * e.g. C(1,2,:) shows activity in neuron 2 following an event in neuron 1
 * C is [Neurons x N eurons x time] with 0 as the first time lag
 * To get the pearson correlation use: squeeze(C(2,1,:)./sqrt(C(1,1,1).*C(2,2,1))));
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

/********************************************************************************************/
/*This function defines the incoming and outgoing matlab variables
 * and sets pointers to the begining of each one. Then calls to the
 * function that makes the new SBEmatrix and to the correlation calculation
 * (the in comming array must me all in the same unit)*/

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *t, *ic, *maxLag, t1, t2;
    long int tmp, startIdx, neu1, neu2, spk1, spk2, nSpk1, nSpk2;
    long int nSpikes, nRowsIc, nNeurons, nLags, lag;
    mwSize dimsOutput[3];
    double *C;
    
    /*check if format is correct-help*/
    if (nrhs!=3)
    {
        mexPrintf("\nThe number of variables entered is: %d\n",nrhs);
        mexPrintf("The correct execution format is:\n[C]=crossCorrRaster(t,ic,lag);\nt is the time vector\nic is the index channel\nlag is the desired lag\nand C is the calculated cross correlation matrix\ne.g. C(1,2,:) shows activity in neuron 2 following an event in neuron 1\nC is [Neurons x Neurons x time] with 0 as the first time lag\nTo get the pearson correlation use: squeeze(C(2,1,:)./sqrt(C(1,1,1).*C(2,2,1))));\n");
        mexErrMsgTxt("See description above");
    }
    
    /* get input variables and their sizes */
    t = mxGetPr(prhs[0]); /*t - spikes times*/
    ic = mxGetPr(prhs[1]); /*ic (4x number of neurons)*/
    maxLag = mxGetPr(prhs[2]); /*t - spikes times*/
    lag = (long int)maxLag[0];
    
    nSpikes = (long int)mxGetN(prhs[0]);         /* total number of spikes*/
    nRowsIc = (long int)mxGetM(prhs[1]);       /* =4  - return number of rows in array*/
    nNeurons = (long int)mxGetN(prhs[1]);         /* number of neurions - return number of columns in array.*/
    nLags = (long int)mxGetN(prhs[2]);         /* number of neurions - return number of columns in array.*/
    
    if (nRowsIc!=4)
        mexErrMsgTxt("\nThe number of rows in ic (index channel) is not 4\n");
    
    /*Checks that the input varibles are not empty matrices*/
    if (nNeurons==0 || nSpikes==0 || nLags!=1)
        mexErrMsgTxt("\nOne of the input variables is either empty or has the wrong size\n");
    
    /* Create output Matrix and set output pointers :*/
    dimsOutput[0]=(mwSize)nNeurons;
    dimsOutput[1]=(mwSize)nNeurons;
    dimsOutput[2]=(mwSize)lag;
    
    plhs[0] = mxCreateNumericArray(3,dimsOutput,mxDOUBLE_CLASS,mxREAL);
    C = mxGetPr(plhs[0]);
    
    mexPrintf("\nSpike train cross-corr neuron:");
    for (neu1=0;neu1<nNeurons;neu1++)
    {
        mexPrintf("%d,",neu1);
        mexEvalString("drawnow;");
        nSpk1=(long int)(ic[neu1*4+3]-ic[neu1*4+2]+1);
        for (neu2=0;neu2<nNeurons;neu2++)
        {
            nSpk2=(long int)(ic[neu2*4+3]-ic[neu2*4+2]+1);
            startIdx=0;
            for (spk1=startIdx;spk1<nSpk1;spk1++)
                for (spk2=0;spk2<nSpk2;spk2++) /*the spikes are sorted!!! so not all indices have to revisited*/
                {
                    t1=t[(long int)ic[neu1*4+2]-1+spk1];
                    t2=t[(long int)ic[neu2*4+2]-1+spk2];
                    if (t2>=t1)
                        if (t2<(t1+lag))
                        {
                            startIdx=spk1;
                            C[(long int)(t2-t1)*nNeurons*nNeurons + neu2*nNeurons + neu1]++;
                        }
                        else
                            break;
                    
                }
        }
    }
    
} /*void mexFunction*/


