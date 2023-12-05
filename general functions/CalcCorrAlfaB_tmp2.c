/*function to calculate the correlation of a convoluted burst matrix:
[R]=CalcCorrAlfab_tmp2(convM);*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

/* Function calculates the mean value of a string*/

double mean(double str[],int n)
{
    int i;
    double sum=0;
    
    for (i=0;i<n;i++)
    {
        sum+=str[i];
    }
    return (sum/n);
    
}

/* Function calculates the time invariant cross correlation, ans[] must be zeroed*/

int corr(double str1[], double str2[], double ans[], int n)
{
    
    int i,j;
    double var1=0,var2=0,sigma1=0,sigma2=0;
    double mean1,mean2;
    
    mean1=mean(str1,n);
    mean2=mean(str2,n);
    
/* Calculation of the standard variation*/
    
    for (i=0;i<n;i++)
    {
        var1=var1+(str1[i]-mean1)*(str1[i]-mean1);
        var2=var2+(str2[i]-mean2)*(str2[i]-mean2);
    }
    
    sigma1=sqrt(var1);
    sigma2=sqrt(var2);
    
    /*Calculation of the crosscorrelation*/
    
    if (sigma1!=0 && sigma2!=0)
    {
        for (i=0;i<n;i++)
        {
            for (j=0;j<=i;j++)
            {
                ans[i]+=(str1[n+j-i-1]-mean1)*(str2[j]-mean2);
            }
            ans[i]=ans[i]/(sigma1*sigma2);
        }
        
        for (i=0;i<(n-1);i++)
        {
            for (j=0;j<=i;j++)
            {
                ans[2*n-2-i]+=(str2[n+j-i-1]-mean2)*(str1[j]-mean1);
            }
            ans[2*n-2-i]=ans[2*n-2-i]/(sigma1*sigma2);
        }
        return 1;
    }
    else return 0;
}


/*******************************************************************************************************************/
/*This function convlutes two arrays (the response array must be symetric with (2*n+1) cells),
  The output (ans) vector's size is (ndata+nrespns-1)*/

void convlv(double *data, int ndata, double *respns, int nrespns, double *ans)
{
    int i;
    int j;
    
/* Building an array with larger ends, with 'nrespns' zero's on both sides*/
    
    double *dataplus;
    
    for (i=0;i<(nrespns+ndata-1);i++) ans[i]=0;
    
    dataplus = (double*)mxMalloc( sizeof(double)*(ndata + (nrespns*2) - 2) );
    for (i=0;i<(ndata+nrespns*2-2);i++)
    {
        if (i>(nrespns-2)&&i<(ndata+nrespns-1))
        {
            dataplus[i]=data[i-nrespns+1];
        }
        else
        {
            dataplus[i]=0;
        }
    }
    
        /* Calculating convolution*/
    
    for (i=0;i<(nrespns+ndata-1);i++)
    {
        for (j=0;j<nrespns;j++)
        {
            ans[i] = (ans[i])+((respns[j])*(dataplus[i+(nrespns-1)-j]));
        }
        
    }
    mxFree(dataplus);
}



/*************************************************************/
double findmax(double meancorr[],int newidth)
{
    double max=0;
    int i;
    for (i=0;i<newidth;i++)
    {
        if (meancorr[i]>max)
        {
            max=meancorr[i];
        }
    }
    return max;
}




/***************************************************************************************************************/
/*This is the core function. it calculates the time invariat cross corelation function
  for each element of the correlation martix (corrmat)*/

//void CalcCorr (double ***sbematrix, double *corrmat,double *gaussian,
//                           int neunum, int nberst,int width, int newidth, int gaussize)
void CalcCorr (double ***sbematrix, double *corrmat,  int neunum, int nberst,int width)
{
    int i=0,j=0,sbe1=0,sbe2=0,neu=0,activeneu=0,count=0,k=0;
    double maxcorrvalue;
    double *meancorr,*tempcorr,*convec;
    
    meancorr = (double*)mxMalloc(sizeof(double)*(2*width-1));
    tempcorr = (double*)mxMalloc(sizeof(double)*(2*width-1));
    convec = (double*)mxMalloc(sizeof(double)*(width));
    
    
/*Makes the correlation matrix by adding the cross correlations of every two neurons in an SBE
  and fininding the maximum value, the first loop is for the diagonal elements which are always 1*/
    
    for (sbe1=0;sbe1<nberst;sbe1++)
    {
        corrmat[sbe1+(nberst*sbe1)]=1;
    }
    
    
    for (sbe1=1;sbe1<nberst;sbe1++)
    {
        for (sbe2=0;sbe2<sbe1;sbe2++)
        {
/*Zeroing*/
            activeneu=0;
            
            for (i=0;i<(2*width-1);i++)
            {
                meancorr[i]=0;
            }
            
/* Calls for correlation and makes an avarages per neuron between two SBE's*/
            
            for (neu=0;neu<neunum;neu++)
            {
                for (i=0;i<(2*width-1);i++)
                {
                    tempcorr[i]=0;
                }
                
                count=corr(sbematrix[sbe1][neu],sbematrix[sbe2][neu],tempcorr, width);
                
                activeneu=activeneu+count;
                if (count!=0)
                {
                    //										for (i=0;i<(2*newidth-1);i++)
                    for (i=0;i<(2*width-1);i++)
                    {
                        meancorr[i]=meancorr[i]+tempcorr[i];
                    }
                }
                
            }
            
            maxcorrvalue=findmax(meancorr,(width*2-1));
            
            if (activeneu==0)
            {
                maxcorrvalue=0;
                printf("Warning 1 : Correlation (%d,%d) is empty\n",sbe1,sbe2);
            }
            else maxcorrvalue=maxcorrvalue/activeneu;
            
            
            corrmat[sbe1+nberst*sbe2]=maxcorrvalue;
            corrmat[sbe2+nberst*sbe1]=maxcorrvalue;
        }
    }
    mxFree(meancorr);
    mxFree(tempcorr);
    mxFree(convec);
}



/**************************************************************************************/
/*This function defines the incoming and outgoing matlab variables
  and sets pointers to the begining of each one. Then calls to the
  function that claculates the correlation of a convoluted matrix convM
  (the in comming array must me all in the same unit)*/

void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[] )
{
    //  double *index, *times, *locaberst, *gaussian, *sbeM, *ConvSbeM;
    double *ConvSbeM;
    double ***sbematrix;
    
    
    double *corrmat=NULL;             //the returned correlation matrix
    double width;
    int neunum, nberst, i=0,j=0,k=0;
    
/*getting pointers and dimensions of the index matrix and the times
  vector*/
    
    ConvSbeM=mxGetPr(prhs[0]);
    nberst =mxGetDimensions(prhs[0])[0];
    neunum =mxGetDimensions(prhs[0])[1];
    width =mxGetDimensions(prhs[0])[2];
    
    //check if format is correct-help
    if (nrhs!=1) {
        mexPrintf("\nThe number of variables entered is: %d\n", nrhs);
        mexErrMsgTxt("The correct execution format is:\n[R]=CalcCorrAlfab_tmp2(convM);\nconvM is a matrix of size Bind x Neurons x width\n");
    }
    
    /* Set the output pointers :*/
    plhs[0] = mxCreateDoubleMatrix(nberst, nberst, mxREAL);   //return correlation matrix 'corrmat'
    corrmat = mxGetPr(plhs[0]);
    
    sbematrix = (double***)mxMalloc(sizeof(double**)*nberst);
    for (i=0;i<nberst;i++) {
        sbematrix[i] = (double**)mxMalloc(sizeof(double*)*neunum);
        for (j=0;j<neunum;j++)
        {
            
            sbematrix[i][j]=(double*)mxMalloc(sizeof(double)*((int)width));
            for (k=0;k<width;k++)
            {
                sbematrix[i][j][k]=ConvSbeM[i+(j)*nberst+(k)*nberst*neunum];
            }
        }
    }
    
    /*convulate SBE matrix and calculate correlation matrix 'corrmat':*/
    //  CalcCorr (sbematrix, corrmat, gaussian, neunum, nberst, ((int)width), ((int)newidth), gaussize);
    CalcCorr (sbematrix, corrmat, neunum, nberst, ((int)width));
    
    for (i=0;i<nberst;i++)
    {
        for (j=0;j<neunum;j++)
        { 
            mxFree(sbematrix[i][j]);
        }
        mxFree(sbematrix[i]);
    }
    mxFree(sbematrix);
}