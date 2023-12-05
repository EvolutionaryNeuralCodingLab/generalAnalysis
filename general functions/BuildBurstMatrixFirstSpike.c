/*[M]=BuildBurstMatrixFirstSpike(indexChchannel,t,Bind(first:last),width,class);
 * t and Bind must be rounded.
 * a function that returns a matrix of selected bursts.
 * compile function using: mex -largeArrayDims BuildBurstMatrixFirstSpike.c;
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
    double *indexCh, *times, *tTrial;
    double width;
    int matClass;
    mwSize i, rowIdx, nNeurons, nTrials, nSpikes, trial, neu, t, ndims[3];
    
    /*Define pointers to go over array of different classes*/
    uint8_T *Muint8;
    double *Mdouble;
    mxLogical *Mlogical;
    
    /*check if format is correct-help*/
    if ((nrhs<4) | (nrhs>5))
    {
        mexPrintf("\nThe number of variables entered is: %d\n",nrhs);
        mexErrMsgTxt("The correct execution format is:\n[M]=BuildBurstMatrix(indexChchannel,t,startTimes,width,class);\nM is a matrix of size trials x Neurons x width\nBurst window is from [startTimes(i) , startTimes(i)+width]\nWhere all variables are ROUNDED and given in the same units\nClass is the variable class (optional): 1='uint8',2='double',3='logical' (Default 2))");
    }
    
    /* get input variables and their sizes */
    indexCh = mxGetPr(prhs[0]); /*indexChchannel (4x number of neurons)*/
    times = mxGetPr(prhs[1]); /*t - spikes timings*/
    tTrial = mxGetPr(prhs[2]); /*time of beginnings of trials - same units as times*/
    width = (mwSize)mxGetScalar(prhs[3]); /*  width of burst window - same units as times  - totla burst window is [tTrial, tTrial+window]*/
    
    /* get input variables and their sizes */
    rowIdx = (mwSize)mxGetM(prhs[0]);       /* =4  - return number of rows in array*/
    nNeurons = (mwSize)mxGetN(prhs[0]);         /* number of neurions - return number of columns in array.*/
    nSpikes = (mwSize)mxGetN(prhs[1]);         /* total number of spikes*/
    nTrials = (mwSize)mxGetN(prhs[2]);         /* total number of bursts*/
    
    /*check that input variables are integers*/
    for (i=0;i<nSpikes;i++)
    {
        times[i]=(long int)times[i];
    }
    
    for (i=0;i<nTrials;i++)
    {
        tTrial[i]=(long int)tTrial[i];
    }
    
    if (nrhs==5)
    {
        matClass = (int)mxGetScalar(prhs[4]);
    }
    else
    {
        matClass = 2;
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
    if (matClass==1)
    {
        /*mexErrMsgTxt("\nGot out on time\n");*/
        plhs[0] = mxCreateNumericArray(3,ndims,mxUINT8_CLASS,mxREAL);
        Muint8 = mxGetPr(plhs[0]);
        /*mexPrintf("\nInteger class used for output\n");*/
        
        for (neu=0;neu<nNeurons;neu++)
        {
            for (trial=0;trial<nTrials;trial++)
            {
                for (t=((mwSize)indexCh[neu*4+2]-1);t<=((mwSize)indexCh[neu*4+3]-1);t++)
                {
                    
                    if  (   ((times[t])<(tTrial[trial]+width))   &   ((times[t])>=(tTrial[trial]))   )
                    {
                        Muint8[trial+(neu)*nTrials+(((mwSize)times[t])-((mwSize)tTrial[trial]))*nTrials*nNeurons]++;
                        break;
                    }
                    
                }
            }
        }
    }
    else if (matClass==2)
    {
        plhs[0] = mxCreateNumericArray(3,ndims,mxDOUBLE_CLASS,mxREAL);
        Mdouble = mxGetPr(plhs[0]);
        /*mexPrintf("\nDouble class used for output\n");*/
        
        for (neu=0;neu<nNeurons;neu++)
        {
            for (trial=0;trial<nTrials;trial++)
            {
                for (t=((mwSize)indexCh[neu*4+2]-1);t<=((mwSize)indexCh[neu*4+3]-1);t++)
                {
                    
                    if  (   ((times[t])<(tTrial[trial]+width))   &   ((times[t])>=(tTrial[trial]))   )
                    {
                        Mdouble[trial+(neu)*nTrials+(((mwSize)times[t])-((mwSize)tTrial[trial]))*nTrials*nNeurons]++;
                        break;
                    }
                }
            }
        }
        
    }
    else if (matClass==3)
    {
        plhs[0] = mxCreateLogicalArray(3,ndims);
        Mlogical = mxGetPr(plhs[0]);
        
        for (neu=0;neu<nNeurons;neu++)
        {
            for (trial=0;trial<nTrials;trial++)
            {
                for (t=((mwSize)indexCh[neu*4+2]-1);t<=((mwSize)indexCh[neu*4+3]-1);t++)
                {
                    
                    if  (   ((times[t])<(tTrial[trial]+width))   &   ((times[t])>=(tTrial[trial]))   )
                    {
                        Mlogical[trial+(neu)*nTrials+(((mwSize)times[t])-((mwSize)tTrial[trial]))*nTrials*nNeurons]=1;
                        break;
                    }
                }
            }
        }
    }
    
    
    /*for (i=0;i<(nTrials-1+(nNeurons-1)*nTrials+(width-1)*nTrials*nNeurons);i++)
     * M[i]=0;*/
    
}
