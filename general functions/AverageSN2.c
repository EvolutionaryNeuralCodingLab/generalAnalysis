/*[Avg,Std,IsNoise]=AverageSN(Signal,N,Nb,C)
Calculates a moving average
N must be an odd number
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"


/*a function that Calculates the moving average*/
void CalcMovingAverage(double *Signal, double *N, double *Nb, double *C, double *IsNoise, double *Noise, double *STDNoise, long int length)
{
    //Definition and initiation of variables
    int i=0;
    int j=0;
    int index=0;
    int index2=1;
    int mBufferSize=(int)N[0]+1; //Mistake in analyzer-For historical reasons N in the analyzer is defined one less than normal
    int mExtendedBufferSize=(int)Nb[0];
    int mOutLinesCounter=0;
    double avg=0.0;
    double mAverage=0.0;
    double mPrevStdDev=0.0;
    double mStdDiv=0.0;
    double mSum=0.0;
    double mSqrSum=0.0;
    double mLimit=0.0;
    double mPrevLimit=0.0;
    double num;
    double *mBuffer;
    bool *mBufferNoise;
    bool isCurrentRowNoise;
    
    //Allocations in memory
    mBuffer = mxMalloc(mBufferSize*sizeof(double));
    mBufferNoise = mxMalloc(mBufferSize*sizeof(bool));
    
    //Intitiaion of arrays
    for (j=0;j<mBufferSize;j++)
        mBufferNoise[j]=true;
    
    //Main loop
    for (i=0;i<length;i++)
    {
        if(i<mBufferSize-1)
        {
            if(i==0)
            {
                mAverage=Signal[i];
                mPrevStdDev = 0.0;
                mBuffer[i]=Signal[i];
                mSum=Signal[i];
                mSqrSum= pow(mBuffer[i],2);
                //Noise[i]=mAverage;
                //STDNoise[i]=0;
                IsNoise[i]=0; //Since there is an inherent deviation in the first value of the standard deviation
                continue;
            }
            mBuffer[i]=Signal[i];
            avg=mAverage;
            mStdDiv = mPrevStdDev;
            mSum+= Signal[i];
            mSqrSum+= pow(mBuffer[i],2);
            mAverage=mSum/(i+1);
             
            mPrevStdDev=sqrt(mSqrSum/(i+1)-pow(mSum/(i+1),2));
            mLimit=mAverage+3*mPrevStdDev+C[0];
            mPrevLimit=mLimit;
            
            mOutLinesCounter=0;
            if ((i-(mBufferSize-1)/2)>=0)
            {
                Noise[i-(mBufferSize-1)/2]=avg;
                STDNoise[i-(mBufferSize-1)/2]=mPrevStdDev;
            }
            IsNoise[i]=1;
            continue;
        }
        else
        {
            avg = mAverage;
            
            index=i%(mBufferSize);
            index2=(i+1)%(mBufferSize);
            
            mStdDiv = mPrevStdDev;
            
            isCurrentRowNoise = Signal[i]<mPrevLimit;
            mBufferNoise[index]=isCurrentRowNoise;
            mBuffer[index]=Signal[i];
            
            mSum+= mBufferNoise[index] ? mBuffer[index] : 0.0;
            mSum-= mBufferNoise[index2] ? mBuffer[index2] : 0.0;
            mSqrSum+= mBufferNoise[index] ? pow(mBuffer[index],2) : 0.0;
            mSqrSum-= mBufferNoise[index2] ? pow(mBuffer[index2],2) : 0.0;
            mOutLinesCounter+= mBufferNoise[index] ? 0 : 1;
            mOutLinesCounter-= mBufferNoise[index2] ? 0 : 1;
            
            mAverage = mSum/(mBufferSize-1-mOutLinesCounter);
            mPrevStdDev=sqrt((mSqrSum/(mBufferSize-1-mOutLinesCounter))-pow(mAverage,2));
            mPrevLimit=mAverage+3*mPrevStdDev+C[0];
            
            if(!mBufferNoise[index] && mOutLinesCounter>mExtendedBufferSize)
            {
                //mexPrintf("\nReset in line: %d",i);
                mOutLinesCounter=0;   
                mSqrSum = 0;
                mSum = 0;
                
                mBufferNoise[(i - 0 + 1)%(mBufferSize)]=true; //To compensate the for reset noise raws
                
                for(j=1;j<mBufferSize;j++)
                {
                    mBufferNoise[(i - j + 1)%(mBufferSize)]=true;
                    num = mBuffer[(i - j + 1)%(mBufferSize)];
                    mSqrSum+=pow(num,2);
                    mSum+=num;
                }
                
                mAverage = mSum/(mBufferSize-1- mOutLinesCounter);
                mPrevStdDev=sqrt((mSqrSum/(mBufferSize-1-mOutLinesCounter))-pow((mSum/(mBufferSize-1-mOutLinesCounter)),2));
                mLimit=mAverage+3*mPrevStdDev+C[0];
                mPrevLimit = mLimit;
            }
            IsNoise[i]=isCurrentRowNoise;
            Noise[i-(mBufferSize-1)/2]=avg;
            STDNoise[i-(mBufferSize-1)/2]=mPrevStdDev;
            
            if (i>(length-mBufferSize)) //temporal patch to fix the right boundry of the signal
            {
                Noise[i]=avg;
                STDNoise[i]=mPrevStdDev;
            }
            continue;
        }
    }
    mxFree(mBuffer);
    mxFree(mBufferNoise);
}
    



/********************************************************************************************/
/*This function defines the incoming and outgoing matlab variables
  and sets pointers to the begining of each one. Then calls to the
  function that makes the new SBEmatrix and to the correlation calculation
  (the in comming array must me all in the same unit)*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
   double *N, *Nb, *C, *Signal;
   long int length=0;
   double *STDNoise=NULL,*Noise=NULL,*IsNoise=NULL;  
    
/*getting pointers and dimensions of the index matrix and the times 
  vector*/

  if (nrhs!=4)
  {
      mexPrintf("\nThe number of variables entered is: %d\n",nrhs);
      mexErrMsgTxt("The correct format is:\n[Avg,Std]=AverageSN(Signal,N,Nb,C)\n");
  }   
  
  Signal = mxGetPr(prhs[0]);            //Signal
  N = mxGetPr(prhs[1]);                 //Size of average window 
  Nb = mxGetPr(prhs[2]);                //Number of allowed unnoisy signals in window
  C = mxGetPr(prhs[3]);                 //Bias
  length = (long int)mxGetDimensions(prhs[0])[1];  //Signal's length
  if (length==1)
      length = (long int)mxGetDimensions(prhs[0])[0];

  /* Set the output pointers :*/
  plhs[0] = mxCreateDoubleMatrix(1, length, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, length, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, length, mxREAL);
  Noise=mxGetPr(plhs[0]);
  STDNoise=mxGetPr(plhs[1]);
  IsNoise=mxGetPr(plhs[2]);
  
  CalcMovingAverage(Signal, N, Nb, C, IsNoise, Noise, STDNoise,length);
  //mexPrintf("\nFinished Successfully\n");
}
