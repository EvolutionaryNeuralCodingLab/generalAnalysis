/*[M]=BuildBurstMatrix(indexchannel,t,Bind(first:last),width);
 t and Bind must be rounded.
a function that returns a matrix of selected bursts.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"


/*a function that builds the SBE matrix*/
void MakeSBE (double *index, double *sbematrix, double *times, double *locaburst,
int neunum, long int ntimes, int nburst, double width)
{
    int t=0,currentburst,neu=0;
    for (neu=0;neu<neunum;neu++)    //run along each neuron:
    {
        //        printf('first: %d \n',(int)index[neu*4+2]-1);
        for (t=((long int)index[neu*4+2]-1);t<=((long int)index[neu*4+3]-1);t++) //t=indexchannel(3,i):indexxchannel(4,i)
        {
            //for each spike, find if it is inside a burst:
            for (currentburst=0;currentburst<nburst;currentburst++)
            {
                //check if: tBurst<t<tBurst+window
                if  (   ((times[t])<(locaburst[currentburst]+width))   &   ((times[t])>=(locaburst[currentburst]))   )
                    sbematrix[currentburst+(neu)*nburst+(((long int)times[t])-((long int)locaburst[currentburst]))*nburst*neunum]++;
            }
        }
    }
}



/********************************************************************************************/
/*This function defines the incoming and outgoing matlab variables
  and sets pointers to the begining of each one. Then calls to the
  function that makes the new SBEmatrix and to the correlation calculation
  (the in comming array must me all in the same unit)*/

void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[] )
{
    double *index, *times, *locaburst, *sbeM;
    double width;     //burst window is from [time , time+width]
    long int ntimes;
    
    int ndims[3], rowindex, neunum, nburst, i=0;
    
/*the width is definded as half of the SBE*/
    
/*getting pointers and dimensions of the index matrix and the times
  vector*/
    
//check if format is correct-help
    if (nrhs!=4)
    {
        mexPrintf("\nThe number of variables entered is: %d\n",nrhs);
        mexErrMsgTxt("The correct execution format is:\n[M]=BuildBurstMatrix(indexchannel,t,Bind,width);\nM is a matrix of size Bind x Neurons x width\nBurst window is from [Bind(i) , Bind(i)+width]\nWhere all variables are ROUNDED and given in the same units\n");
    }
    
    index = mxGetPr(prhs[0]); //indexchannel (4x number of neurons)
    times = mxGetPr(prhs[1]); //  t - spikes timings
    locaburst = mxGetPr(prhs[2]); //location of beginning of burst window - same units as times
    width = mxGetScalar(prhs[3]); //  width of burst window - same units as times  - totla burst window is [locaburst, locaburst+window]

    
    rowindex = mxGetM(prhs[0]);       // =4  - return number of rows in array
    neunum = mxGetN(prhs[0]);         // number of neurions - return number of columns in array.
    ntimes = mxGetN(prhs[1]);         // total number of spikes
    nburst = mxGetN(prhs[2]);         // total number of bursts

    //Checks that the input varibles are not empty matrices
    if (neunum==0 | ntimes==0 | nburst==0)
    {
        mexErrMsgTxt("\nOne of the input variables is empty\n");
    }
    
/*After convolution the vector size is bigger by the size of (the
  response function(respns))-1)*/
    
/*defenition of a 3D matrix in dynamic memory*/
    
  /*dimentions of burst matrixes: sbeM, ConvSbeM:*/
    ndims[0]=nburst;
    ndims[1]=neunum;
    ndims[2]=(int)width;
    
  /* Set the output pointers :*/
    plhs[0] = mxCreateNumericArray(3,ndims,mxDOUBLE_CLASS,mxREAL);  //return burst matrix (sbeM)
    sbeM=mxGetPr(plhs[0]);
    for (i=0;i<(nburst-1+(neunum-1)*nburst+(width-1)*nburst*neunum);i++)
        sbeM[i]=0;
    
  /* Error Messages */
    if (locaburst[0]<0) mexWarnMsgTxt("Burst location starts from a negative value\n");
    
    //build SBE matrix:
    MakeSBE (index, sbeM, times, locaburst, neunum, ntimes, nburst, width);
}
