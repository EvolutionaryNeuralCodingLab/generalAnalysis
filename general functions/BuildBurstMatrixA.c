/*[M]=BuildBurstMatrix(indexChchannel,t,Bind(first:last),width);
 * t and Bind must be rounded.
 * a function that returns a matrix of selected bursts.
 * compile function using: mex -largeArrayDims BuildBurstMatrix.c;
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
    double *indexCh, *times, *tTrial, *intensities;
    double width;
    mwSize i, rowIdx, nNeurons, nTrials, nSpikes, trial, neu, t, ndims[3];
    
    /*Define pointers to go over array of different classes*/
    uint8_T *Muint8;
    double *Mdouble;
    mxLogical *Mlogical;
    
    /*check if format is correct-help*/
    if (nrhs!=5)
    {
        mexPrintf("\nThe number of variables entered is: %d\n",nrhs);
        mexErrMsgTxt("The correct execution format is:\n[M]=BuildBurstMatrixA(ic,t,I,Bind,width,class);\nM is a  matrix of size Bind x Neurons x width\nBurst window is from [Bind(i) , Bind(i)+width]\nWhere all variables are ROUNDED and given in the same units\n");
    }
    
    /* get input variables and their sizes */
    indexCh = mxGetPr(prhs[0]); /*indexChchannel (4x number of neurons)*/
    times = mxGetPr(prhs[1]); /*t - spikes timings*/
    intensities = mxGetPr(prhs[2]); /*  I - the intensities in times t*/
    tTrial = mxGetPr(prhs[3]); /*time of beginnings of trials - same units as times*/
    width = (mwSize)mxGetScalar(prhs[4]); /*  width of burst window - same units as times  - totla burst window is [tTrial, tTrial+window]*/
    
    
    rowIdx = (mwSize)mxGetM(prhs[0]);       /* =4  - return number of rows in array*/
    nNeurons = (mwSize)mxGetN(prhs[0]);         /* number of neurions - return number of columns in array.*/
    nSpikes = (mwSize)mxGetN(prhs[1]);         /* total number of spikes*/
    nTrials = (mwSize)mxGetN(prhs[3]);         /* total number of bursts*/
    
    /*check that input variables are integers*/
    for (i=0;i<nSpikes;i++)
    {
        times[i]=(long int)times[i];
    }
    
    for (i=0;i<nTrials;i++)
    {
        tTrial[i]=(long int)tTrial[i];
    }
    

    /*Checks that the input varibles are not empty matrices*/
    if (nNeurons==0 || nSpikes==0 || nTrials==0)
    {
        mexErrMsgTxt("\nOne of the input variables is empty\n");
    }
    
    /*After convolution the vector size is bigger by the size of (the
     * response function(respns))-1)*/
    
    /*dimentions of burst matrixes: M, ConvSbeM:*/
    ndims[0]=nTrials;
    ndims[1]=nNeurons;
    ndims[2]=(mwSize)width;
    
    /* Error Messages */
    if (tTrial[0]<0)
        mexWarnMsgTxt("Burst location starts from a negative value\n");
    
    /* Set the output pointers :*/
 

        plhs[0] = mxCreateNumericArray(3,ndims,mxDOUBLE_CLASS,mxREAL);
        Mdouble = mxGetPr(plhs[0]);
        /*mexPrintf("\nDouble class used for output\n");*/
        
        for (neu=0;neu<nNeurons;neu++)
        {
            for (t=((mwSize)indexCh[neu*4+2]-1);t<=((mwSize)indexCh[neu*4+3]-1);t++)
            {
                for (trial=0;trial<nTrials;trial++)
                {
                    if  (   ((times[t])<(tTrial[trial]+width))   &   ((times[t])>=(tTrial[trial]))   )
                        Mdouble[trial+(neu)*nTrials+(((mwSize)times[t])-((mwSize)tTrial[trial]))*nTrials*nNeurons] += intensities[t];
                }
            }
        }
        
}
