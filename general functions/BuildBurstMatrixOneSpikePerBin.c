/*[M]=BuildBurstMatrix(indexchannel,t,Bind(first:last),width);
a function that returns a matrix of selected bursts.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"


/*a function that builds the SBE matrix*/
void MakeSBE (double *index, double ***sbematrix, double *times, double *locaberst, 
              int neunum, long int ntimes, int nberst, double width)
{
        int t=0,currentberst,neu=0;
        for (neu=0;neu<neunum;neu++)    //run along each neuron:
        {
//        printf('first: %d \n',(int)index[neu*4+2]-1);
				for (t=((long int)index[neu*4+2]-1);t<=((long int)index[neu*4+3]-1);t++) //t=indexchannel(3,i):indexxchannel(4,i)
                {
                        //for each spike, find if it is inside a burst:
                        for (currentberst=0;currentberst<nberst;currentberst++)        
                         {
                                //check if: tBurst<t<tBurst+window 
                                if (  ((times[t])<(locaberst[currentberst]+width)) &  ((times[t])>=(locaberst[currentberst]))   ) 
                                                //add +1 in SBE matrix:
                                                sbematrix[currentberst][neu][((long int)times[t])-((long int)locaberst[currentberst])]=1;
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
   double *index, *times, *locaberst, *sbeM;
   double ***sbematrix;
   
 double *corrmat=NULL;
  double width;     //burst window is from [time , time+width]
  long int ntimes;

    int rowindex, neunum, nberst, i=0, j=0, k=0;    
    
int ndims[3];
  
/*the width is definded as half of the SBE*/  
  
/*getting pointers and dimensions of the index matrix and the times 
  vector*/
  
  index = mxGetPr(prhs[0]); //indexchannel (4x number of neurons)
  times = mxGetPr(prhs[1]); //  t - spikes timings
  locaberst = mxGetPr(prhs[2]); //location of beginning of burst window - same units as times
  width = mxGetScalar(prhs[3]); //  width of burst window - same units as times  - totla burst window is [locaburst, locaburst+window]
  
  rowindex = mxGetM(prhs[0]);       // =4  - return number of rows in array
  neunum = mxGetN(prhs[0]);         // number of neurions - return number of columns in array.
  ntimes = mxGetN(prhs[1]);         // total number of spikes
  nberst = mxGetN(prhs[2]);         // total number of bursts
  
/*After convolution the vector size is bigger by the size of (the 
  response function(respns))-1)*/

/*defenition of a 3D matrix in dynamic memory*/

//put zeros in SBE matrix:
  sbematrix = (double***)mxMalloc(sizeof(double**)*nberst);
  for (i=0;i<nberst;i++)
  {
          sbematrix[i] = (double**)mxMalloc(sizeof(double*)*neunum);
          for (j=0;j<neunum;j++)
          {
                    sbematrix[i][j]=(double*)mxMalloc(sizeof(double)*((int)width));
                  for (k=0;k<width;k++) 
                  {
                          sbematrix[i][j][k]=0;
                  }
          }
  }
  
  /*dimentions of burst matrixes: sbeM, ConvSbeM:*/
  ndims[0]=nberst;
  ndims[1]=neunum;
  ndims[2]=(int)width;

  /* Set the output pointers :*/
    plhs[0] = mxCreateNumericArray(3,ndims,mxDOUBLE_CLASS,mxREAL);  //return burst matrix (sbeM)
  sbeM=mxGetPr(plhs[0]);


  /* Error Messages */
//  if (locaberst[nberst-1]>times[ntimes-1]) printf("Error 1 : berst location exceeds recorded data\n");
  if (locaberst[0]<0) printf("Error 2 : Berst location starts from a negative value\n");

//build SBE matrix:  
  MakeSBE (index, sbematrix, times, locaberst, neunum, ntimes, nberst, width);
  // copy burst matrix to returned pointer:
  for (i=0;i<nberst;i++)
    {
           for (j=0;j<neunum;j++)
           {
                  for (k=0;k<width;k++) 
                  {
                  //matlab accepts matrixes as vectors:
                          sbeM[i+(j)*nberst+(k)*nberst*neunum]=sbematrix[i][j][k];
                  }
                  mxFree(sbematrix[i][j]);
          }
          mxFree(sbematrix[i]);
    }
    mxFree(sbematrix); 
}
