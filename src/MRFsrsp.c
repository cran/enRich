//***********************************************************//
//***********************************************************//
//***    The program is used to do ChIP-Seq data analysis ***//
//***    for histone modification using Markov Random     ***//
//***    Field (MRF model), which will form the core part ***//
//***    of package RFseq. The program is written by      ***//
//***    Dr Yanchun Bao, of Brunel university.            ***//
//***********************************************************//
//***********************************************************//

//#include <iostream.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#define cp 5 // the cutting point for determine the initial states  
#define I 2 // dimension of state space of X


//**********************************************************//
//*********************   Main program  ********************//
void MRFsrsp(int *data1, int *data2, int *nsr, int *nsp, int *size, int *met, double *qprior1, double *piprior1, double *poisprior1, double *NBprior1, double *qprior2, double *piprior2, double *poisprior2, double *NBprior2, int *N, int *Nb, int *jumpN, double *var, double *var1, double *var2, double *PP1, double *PP2, double *es_q11, double *es_q10, double *es_lambda11, double *es_mu11, double *es_phi11, double *es_pi1, double *es_lambda10, double *es_mu10, double *es_phi10, double *es_q21, double *es_q20, double *es_lambda21, double *es_mu21, double *es_phi21, double *es_pi2, double *es_lambda20, double *es_mu20, double *es_phi20, double *loglikeli, double *acrate, double *acrate1, double *acrate2)
{
  int S=*size;
  int N1=*N;
  int bN=*Nb;
  int jN=*jumpN;
  int nsr1=*(nsr+0);
  int nsr2=*(nsr+1);
  int nsp1=*(nsp+0);
  int nsp2=*(nsp+1);
  int nr1=nsr1+nsp1;
  int nr2=nsr2+nsp2;
  int nr1I=nr1*I;
  int nr2I=nr2*I;
  int nr1S=nr1*S;
  int nr2S=nr2*S;
  int nsp1I=nsp1*I;
  int nsp2I=nsp2*I;
  int nsp1S=nsp1*S;
  int nsp2S=nsp2*S;
  //Rprintf("S, N and nsr, nsp are: %d %d %d %d %d %d %d %d\n", S, N1, nsr1, nsr2, nsp1, nsp2, nr1, nr2);
  int n[I][I];//the counts of pair
  double Q[I][I];
  int i, i1, i2, j, j1, j2, s, z; 
  double *sum2X;//using for counting pairs of X
  sum2X=(double*)malloc(sizeof(double*)*(S-1));
  double U,U1,temp, temp1, temp2,a,b,logl;

  // Definition used for first data set of condition 1
  double AQrep1, BQrep1, *AQnrep1, *BQnrep1;  // priors of transition probabilities q0.
  double *Api1, *Bpi1;// priors of mixture portion of ZIP model.
  Api1=(double*)malloc(sizeof(double*)*nr1);
  Bpi1=(double*)malloc(sizeof(double*)*nr1);
  double *Alambda1, *Blambda1; // priors of lambda for poisson fitting.
  Alambda1=(double*)malloc(sizeof(double*)*nr1I);
  Blambda1=(double*)malloc(sizeof(double*)*nr1I);
  double *Amu1, *Bmu1;//priors of mu_i, which is the mean of NB for ith component;
  Amu1=(double*)malloc(sizeof(double*)*nr1I);
  Bmu1=(double*)malloc(sizeof(double*)*nr1I);
  double *Aphi1, *Bphi1;//priors of phi_i, which is the overdispersion of NB for ithe component;
  Aphi1=(double*)malloc(sizeof(double*)*nr1I);
  Bphi1=(double*)malloc(sizeof(double*)*nr1I);
  int *dataset1;// bin counts--observations.
  dataset1=(int*)malloc(sizeof(int*)*nr1S);
  int *stateXrep1; // bin states, latent variables for replicates
  int *stateXnrep1; // bin states, latent variables for nonreplicates
  int *stateZ1; // inner latent variables
  stateZ1=(int*)malloc(sizeof(int*)*nr1S);
  double *lambda1, *mu1, *phi1, *pi1, *q1, *qrep1, *qnrep1; // parameters of Poisson-Gamma mixing distribution and transition probabilities
  lambda1=(double*)malloc(sizeof(double*)*nr1I);
  mu1=(double*)malloc(sizeof(double*)*nr1I);
  phi1=(double*)malloc(sizeof(double*)*nr1I);
  pi1=(double*)malloc(sizeof(double*)*nr1);
  q1=(double*)malloc(sizeof(double*)*nr1I);
  double *probXrep1, *sumprobXrep1, *probXnrep1, *sumprobXnrep1, *probZ1;
  int *indexY1;
  probZ1=(double*)malloc(sizeof(double*)*nr1); 
  indexY1=(int*)malloc(sizeof(int*)*nr1S);
  double **logpY1,**pY1,**pXcYrep1, **pXcYnrep1; 
  logpY1=(double**)malloc(sizeof(double*)*I);
  pY1=(double**)malloc(sizeof(double*)*I);
  for(i=0;i<I;i++)
    {
      logpY1[i]=(double*)malloc(sizeof(double*)*nr1S);
      pY1[i]=(double*)malloc(sizeof(double*)*nr1S);
    }
  int indexnsr1;
  if (nsr1>0)
    {
      indexnsr1=1;
      stateXrep1=(int*)malloc(sizeof(int*)*S);//********************************check********************************************************
      qrep1=(double*)malloc(sizeof(double*)*I);
      probXrep1=(double*)malloc(sizeof(double*)*S);//********************************check********************************************************
      sumprobXrep1=(double*)malloc(sizeof(double*)*S);//********************************check********************************************************
      pXcYrep1=(double**)malloc(sizeof(double*)*I);
      for (i=0;i<I;i++)
	{
	  pXcYrep1[i]=(double*)malloc(sizeof(double*)*S);//********************************check********************************************************
	}      
    }
  else
    {
      indexnsr1=0;
      AQrep1=0;
      BQrep1=0;
      stateXrep1=NULL;
      qrep1=NULL;
      probXrep1=NULL;
      sumprobXrep1=NULL;
      pXcYrep1=NULL;
    }
  if (nsp1>0)
    {
      AQnrep1=(double*)malloc(sizeof(double*)*nsp1);
      BQnrep1=(double*)malloc(sizeof(double*)*nsp1);
      stateXnrep1=(int*)malloc(sizeof(int*)*nsp1S);//********************************check********************************************************
      qnrep1=(double*)malloc(sizeof(double*)*nsp1I);//********************************check********************************************************
      probXnrep1=(double*)malloc(sizeof(double*)*nsp1S);//********************************check********************************************************
      sumprobXnrep1=(double*)malloc(sizeof(double*)*nsp1S);//********************************check********************************************************
      pXcYnrep1=(double**)malloc(sizeof(double*)*I);
      for(i=0;i<I;i++)
	{
	  pXcYnrep1[i]=(double*)malloc(sizeof(double*)*nsp1S);//********************************check********************************************************
	}
    }
  else
    {
      AQnrep1=NULL;
      BQnrep1=NULL;
      stateXnrep1=NULL;
      qnrep1=NULL;
      probXnrep1=NULL;
      sumprobXnrep1=NULL;
      pXcYnrep1=NULL;
    }

  // Definition used for second data set of condition 2
  double AQrep2, BQrep2, *AQnrep2, *BQnrep2;  // priors of transition probabilities q0.
  double *Api2, *Bpi2;// priors of mixture portion of ZIP model.
  Api2=(double*)malloc(sizeof(double*)*nr2);
  Bpi2=(double*)malloc(sizeof(double*)*nr2);
  double *Alambda2, *Blambda2; // priors of lambda for poisson fitting.
  Alambda2=(double*)malloc(sizeof(double*)*nr2I);
  Blambda2=(double*)malloc(sizeof(double*)*nr2I);
  double *Amu2, *Bmu2;//priors of mu_i, which is the mean of NB for ith component;
  Amu2=(double*)malloc(sizeof(double*)*nr2I);
  Bmu2=(double*)malloc(sizeof(double*)*nr2I);
  double *Aphi2, *Bphi2;//priors of phi_i, which is the overdispersion of NB for ithe component;
  Aphi2=(double*)malloc(sizeof(double*)*nr2I);
  Bphi2=(double*)malloc(sizeof(double*)*nr2I);
  int *dataset2;// bin counts--observations.
  dataset2=(int*)malloc(sizeof(int*)*nr2S);
  int *stateXrep2; // bin states, latent variables for replicates
  int *stateXnrep2; // bin states, latent variables for nonreplicates
  int *stateZ2; // inner latent variables
  stateZ2=(int*)malloc(sizeof(int*)*nr2S);
  double *lambda2, *mu2, *phi2, *pi2, *q2, *qrep2, *qnrep2; // parameters of Poisson-Gamma mixing distribution and transition probabilities
  lambda2=(double*)malloc(sizeof(double*)*nr2I);
  mu2=(double*)malloc(sizeof(double*)*nr2I);
  phi2=(double*)malloc(sizeof(double*)*nr2I);
  pi2=(double*)malloc(sizeof(double*)*nr2);
  q2=(double*)malloc(sizeof(double*)*nr2I);
  double *probXrep2, *sumprobXrep2, *probXnrep2, *sumprobXnrep2, *probZ2;
  int *indexY2;
  probZ2=(double*)malloc(sizeof(double*)*nr2); 
  indexY2=(int*)malloc(sizeof(int*)*nr2S);
  double **logpY2,**pY2,**pXcYrep2, **pXcYnrep2; 
  logpY2=(double**)malloc(sizeof(double*)*I);
  pY2=(double**)malloc(sizeof(double*)*I);
  for(i=0;i<I;i++)
    {
      logpY2[i]=(double*)malloc(sizeof(double*)*nr2S);
      pY2[i]=(double*)malloc(sizeof(double*)*nr2S);
    }
  int indexnsr2;
  if (nsr2>0)
    {
      indexnsr2=1;
      stateXrep2=(int*)malloc(sizeof(int*)*S);//********************************check********************************************************
      qrep2=(double*)malloc(sizeof(double*)*I);
      probXrep2=(double*)malloc(sizeof(double*)*S);//********************************check********************************************************
      sumprobXrep2=(double*)malloc(sizeof(double*)*S);//********************************check********************************************************
      pXcYrep2=(double**)malloc(sizeof(double*)*I);
      for (i=0;i<I;i++)
	{
	  pXcYrep2[i]=(double*)malloc(sizeof(double*)*S);//********************************check********************************************************
	}      
    }
  else
    {
      indexnsr2=0;
      AQrep2=0;
      BQrep2=0;
      stateXrep2=NULL;
      qrep2=NULL;
      probXrep2=NULL;
      sumprobXrep2=NULL;
      pXcYrep2=NULL;
    }
  if (nsp2>0)
    {
      AQnrep2=(double*)malloc(sizeof(double*)*nsp2);
      BQnrep2=(double*)malloc(sizeof(double*)*nsp2);
      stateXnrep2=(int*)malloc(sizeof(int*)*nsp2S);//********************************check********************************************************
      qnrep2=(double*)malloc(sizeof(double*)*nsp2I);//********************************check********************************************************
      probXnrep2=(double*)malloc(sizeof(double*)*nsp2S);//********************************check********************************************************
      sumprobXnrep2=(double*)malloc(sizeof(double*)*nsp2S);//********************************check********************************************************
      pXcYnrep2=(double**)malloc(sizeof(double*)*I);
      for(i=0;i<I;i++)
	{
	  pXcYnrep2[i]=(double*)malloc(sizeof(double*)*nsp2S);//********************************check********************************************************
	}
    }
  else
    {
      AQnrep2=NULL;
      BQnrep2=NULL;
      stateXnrep2=NULL;
      qnrep2=NULL;
      probXnrep2=NULL;
      sumprobXnrep2=NULL;
      pXcYnrep2=NULL;
    }

  double *pX1, *pX0;
  int nsrsp=indexnsr1+indexnsr2+nsp1+nsp2;
  pX1=(double*)malloc(sizeof(double*)*nsrsp);
  pX0=(double*)malloc(sizeof(double*)*nsrsp);
   
  ///////// R data1 pass to data
  for(s=0;s<nr1S;s++)
    {
      dataset1[s]=*(data1+s);     
    }
   for(s=0;s<nr2S;s++)
    {
      dataset2[s]=*(data2+s); 
    }

/*
   for(s=15;s<20;s++)
   {
     Rprintf("dataset1[:,]= %d %d\n", dataset1[s], dataset1[(nr1-1)*S+s]);
   }
   for(s=15;s<20;s++)
   {
     Rprintf("dataset2[:,]= %d %d\n", dataset2[s], dataset2[(nr1-1)*S+s]);
   }
*/

  for (s=0;s<S;s++)
    {
      if (nsr1>0)
	{
	  sumprobXrep1[s]=0.0;
	}
      if (nsr2>0)
	{
	  sumprobXrep2[s]=0.0;
	}
    }
  if (nsp1>0)
    {
      for (s=0;s<nsp1S;s++)
	{
	  sumprobXnrep1[s]=0.0;
	}
    }
  if (nsp2>0)
    {
      for (s=0;s<nsp2S;s++)
	{
	  sumprobXnrep2[s]=0.0;
	}
    }

  ///////// give priors for parameters; 
  if(nsr1>0)
    {
      AQrep1=*(qprior1+0);//qprior1
      BQrep1=*(qprior1+1); 
    }
  if(nsp1>0)
    {
      for (i=0;i<nsp1;i++)
	{
	  AQnrep1[i]=*(qprior1+i*I+I*indexnsr1+0);
	  BQnrep1[i]=*(qprior1+i*I+I*indexnsr1+1);
	} 
    }    
  if (nsr2>0)
    {
      AQrep2=*(qprior2+0);
      BQrep2=*(qprior2+1);
    }
  if (nsp2>0)
    {
      for (i=0;i<nsp2;i++)
	{
	  AQnrep2[i]=*(qprior2+i*I+I*indexnsr2+0);
	  BQnrep2[i]=*(qprior2+i*I+I*indexnsr2+1);
	}
    }
 
  for (i=0;i<nr1;i++)
    {
      Api1[i]=*(piprior1+i*I+0);
      Bpi1[i]=*(piprior1+i*I+1);
  
      Alambda1[i*I+1]=*(poisprior1+i*4+0);
      Blambda1[i*I+1]=*(poisprior1+i*4+1);
      Alambda1[i*I+0]=*(poisprior1+i*4+2);
      Blambda1[i*I+0]=*(poisprior1+i*4+3);

      Amu1[i*I+1]=*(NBprior1+i*8+0);
      Bmu1[i*I+1]=*(NBprior1+i*8+1);
      Aphi1[i*I+1]=*(NBprior1+i*8+2);
      Bphi1[i*I+1]=*(NBprior1+i*8+3);
      Amu1[i*I+0]=*(NBprior1+i*8+4);
      Bmu1[i*I+0]=*(NBprior1+i*8+5);
      Aphi1[i*I+0]=*(NBprior1+i*8+6);
      Bphi1[i*I+0]=*(NBprior1+i*8+7);

      //Rprintf("prior of q= %f %f %f %f \n",AQ[i*I+0], AQ[i*I+1], BQ[i*I+0], BQ[i*I+1]);
      //Rprintf("Prior of %d",i);
      //Rprintf("th experiment, prior of pi= %f %f \n", Api[i], Bpi[i]);
      //Rprintf("prior of lambda= %f %f %f %f \n", Alambda[i*I+1], Blambda[i*I+1], Alambda[i*I+0], Blambda[i*I+0]);
      //Rprintf("prior of mu and phi= %f %f %f %f %f %f %f %f \n", Amu[i*I+1], Bmu[i*I+1], Aphi[i*I+1], Bphi[i*I+1], Amu[i*I+0], Bmu[i*I+0], Aphi[i*I+0], Bphi[i*I+0]); 
    }
  
  for (i=0;i<nr2;i++)
    {
      Api2[i]=*(piprior2+i*I+0);
      Bpi2[i]=*(piprior2+i*I+1);
  
      Alambda2[i*I+1]=*(poisprior2+i*4+0);
      Blambda2[i*I+1]=*(poisprior2+i*4+1);
      Alambda2[i*I+0]=*(poisprior2+i*4+2);
      Blambda2[i*I+0]=*(poisprior2+i*4+3);

      Amu2[i*I+1]=*(NBprior2+i*8+0);
      Bmu2[i*I+1]=*(NBprior2+i*8+1);
      Aphi2[i*I+1]=*(NBprior2+i*8+2);
      Bphi2[i*I+1]=*(NBprior2+i*8+3);
      Amu2[i*I+0]=*(NBprior2+i*8+4);
      Bmu2[i*I+0]=*(NBprior2+i*8+5);
      Aphi2[i*I+0]=*(NBprior2+i*8+6);
      Bphi2[i*I+0]=*(NBprior2+i*8+7);

      //Rprintf("prior of q= %f %f %f %f \n",AQ[i*I+0], AQ[i*I+1], BQ[i*I+0], BQ[i*I+1]);
      //Rprintf("Prior of %d",i);
      //Rprintf("th experiment, prior of pi= %f %f \n", Api[i], Bpi[i]);
      //Rprintf("prior of lambda= %f %f %f %f \n", Alambda[i*I+1], Blambda[i*I+1], Alambda[i*I+0], Blambda[i*I+0]);
      //Rprintf("prior of mu and phi= %f %f %f %f %f %f %f %f \n", Amu[i*I+1], Bmu[i*I+1], Aphi[i*I+1], Bphi[i*I+1], Amu[i*I+0], Bmu[i*I+0], Aphi[i*I+0], Bphi[i*I+0]); 
    }
 
 //////// give initial values, set X=1 if Y>cp, where cp is initial cutting points
  double *sumZrep1, *sumrep1, *sumZnrep1, *sumnrep1;
  int nrep11, nrep10, *nZrep11, *nZrep10, *nnrep11,*nnrep10, *nZnrep11,*nZnrep10;
  if (nsr1>0)
    {
      sumZrep1=(double*)malloc(sizeof(double*)*nsr1);
      sumrep1=(double*)malloc(sizeof(double*)*nsr1);
      nZrep11=(int*)malloc(sizeof(int*)*nsr1);
      nZrep10=(int*)malloc(sizeof(int*)*nsr1);
    }
  else
    {
      sumZrep1=NULL;
      sumrep1=NULL;
      nrep11=0;
      nrep10=0;
      nZrep11=NULL;
      nZrep10=NULL;
    }
  if (nsp1>0)
    {
      sumZnrep1=(double*)malloc(sizeof(double*)*nsp1);
      sumnrep1=(double*)malloc(sizeof(double*)*nsp1);  
      nnrep11=(int*)malloc(sizeof(int*)*nsp1);
      nnrep10=(int*)malloc(sizeof(int*)*nsp1);
      nZnrep11=(int*)malloc(sizeof(int*)*nsp1);
      nZnrep10=(int*)malloc(sizeof(int*)*nsp1);
    }
  else
    {
      sumZnrep1=NULL;
      sumnrep1=NULL;
      nnrep11=NULL;
      nnrep10=NULL;
      nZnrep11=NULL;
      nZnrep10=NULL;
    }
  
  double *sumZrep2, *sumrep2, *sumZnrep2, *sumnrep2;
  int nrep21, nrep20, *nZrep21, *nZrep20, *nnrep21,*nnrep20, *nZnrep21,*nZnrep20;
  if (nsr2>0)
    {
      sumZrep2=(double*)malloc(sizeof(double*)*nsr2);
      sumrep2=(double*)malloc(sizeof(double*)*nsr2);
      nZrep21=(int*)malloc(sizeof(int*)*nsr2);
      nZrep20=(int*)malloc(sizeof(int*)*nsr2);
    }
  else
    {
      sumZrep2=NULL;
      sumrep2=NULL;
      nrep21=0;
      nrep20=0;
      nZrep21=NULL;
      nZrep20=NULL;
    }
  if (nsp2>0)
    {
      sumZnrep2=(double*)malloc(sizeof(double*)*nsp2);
      sumnrep2=(double*)malloc(sizeof(double*)*nsp2);  
      nnrep21=(int*)malloc(sizeof(int*)*nsp2);
      nnrep20=(int*)malloc(sizeof(int*)*nsp2);
      nZnrep21=(int*)malloc(sizeof(int*)*nsp2);
      nZnrep20=(int*)malloc(sizeof(int*)*nsp2);
    }
  else
    {
      sumZnrep2=NULL;
      sumnrep2=NULL;
      nnrep21=NULL;
      nnrep20=NULL;
      nZnrep21=NULL;
      nZnrep20=NULL;
    }
  
  double c=1; //the common portion of transtion probability

  if(nsr1>0)
    {
      nrep11=0;
      nrep10=0;
      for (s=0;s<S;s++)
	{
	  temp1=0;
	  for (i=0;i<nsr1;i++)
	    {
	      temp1=temp1+dataset1[i*S+s];
	    }
	  if (temp1>cp*(nsr1-1))
	    {
	      stateXrep1[s]=1;
	      nrep11=nrep11++;
	    }
	  else
	    {
	      stateXrep1[s]=0;
	      nrep10=nrep10++;
	    }
	}
    }
  if(nsr2>0)
    {
      nrep21=0;
      nrep20=0;
      for (s=0;s<S;s++)
	{
	  temp2=0;
	  for (i=0;i<nsr2;i++)
	    {
	      temp2=temp2+dataset2[i*S+s];
	    }
	  if (temp2>cp*(nsr2-1))
	    {
	      stateXrep2[s]=1;
	      nrep21=nrep21++;
	    }
	  else
	    {
	      stateXrep2[s]=0;
	      nrep20=nrep20++;
	    }
	}
    }
  // set initial values for replicates of first conditon (protein)
  if (nsr1>0)
    {
      for (i=0;i<nsr1;i++)
	{
	  sumZrep1[i]=0.0;
	  sumrep1[i]=0.0;
	  nZrep11[i]=0;
	  nZrep10[i]=0;
	  for (s=0;s<S;s++)
	    { 
	      if(stateXrep1[s]==1)
		{ 
		  sumrep1[i]=sumrep1[i]+dataset1[i*S+s];
		  stateZ1[i*S+s]=2;
		  indexY1[i*S+s]=0;
		}
	      else
		{
		  if (dataset1[i*S+s]==0)
		    {
		      indexY1[i*S+s]=1;
		      stateZ1[i*S+s]=0;
		      nZrep10[i]=nZrep10[i]++;
		    }
		  else
		    {
		      indexY1[i*S+s]=0;
		      stateZ1[i*S+s]=1;
		      sumZrep1[i]=sumZrep1[i]+dataset1[i*S+s];
		      nZrep11[i]=nZrep11[i]++;
		    }
		}
	    }
	  if (sumrep1[i]==0|nrep11==0)
	    {
	      sumrep1[i]=5;
	      nrep11=1;
	    }
	  //printf("initial sumZ1, sum1, sum, n1, n0, n1+n0=% lf %lf %lf %d %d %d \n",sumZ1, sum1, sum, nZ0, nZ1, nZ0+nZ1, n0, n1, n0+n1);
	  if (*met==0)
	    {
	      lambda1[i*I+0]=sumZrep1[i]/nZrep11[i];
	      lambda1[i*I+1]=0;
	      mu1[i*I+1]=sumrep1[i]/nrep11;
	      mu1[i*I+0]=0;
	      phi1[i*I+1]=1;
	      phi1[i*I+0]=0;
	    }
	  if (*met==1)
	    {
	      lambda1[i*I+0]=sumZrep1[i]/nZrep11[i];//sum0/n0;
	      lambda1[i*I+1]=sumrep1[i]/nrep11; 
	      mu1[i*I+0]=0;
	      mu1[i*I+1]=0;
	      phi1[i*I+0]=0;
	      phi1[i*I+1]=0;
	    }
	  if (*met==2)
	    {
	      lambda1[i*I+0]=0;//sum0/n0;
	      lambda1[i*I+1]=0; 
	      mu1[i*I+0]=sumZrep1[i]/nZrep11[i];
	      mu1[i*I+1]=sumrep1[i]/nrep11;
	      phi1[i*I+0]=1;
	      phi1[i*I+1]=1;
	    }
	  pi1[i]=(nZrep11[i]+0.0)/(nZrep10[i]+0.0);
	}
      qrep1[0]=0.1; //prob(X_s=1|X_{s-1}=0)
      qrep1[1]=0.9; //prob(X_s=1|X_{s-1}=1)
    }

  // set initial parameters for replicates of second conditon(protein)
  if (nsr2>0)
    {
      for (i=0;i<nsr2;i++)
	{
	  sumZrep2[i]=0.0;
	  sumrep2[i]=0.0;
	  nZrep21[i]=0;
	  nZrep20[i]=0;
	  for (s=0;s<S;s++)
	    { 
	      if(stateXrep2[s]==1)
		{ 
		  sumrep2[i]=sumrep2[i]+dataset2[i*S+s];
		  stateZ2[i*S+s]=2;
		  indexY2[i*S+s]=0;
		}
	      else
		{
		  if (dataset2[i*S+s]==0)
		    {
		      indexY2[i*S+s]=1;
		      stateZ2[i*S+s]=0;
		      nZrep20[i]=nZrep20[i]++;
		    }
		  else
		    {
		      indexY2[i*S+s]=0;
		      stateZ2[i*S+s]=1;
		      sumZrep2[i]=sumZrep2[i]+dataset2[i*S+s];
		      nZrep21[i]=nZrep21[i]++;
		    }
		}
	    }
	  if (sumrep2[i]==0|nrep21==0)
	    {
	      sumrep2[i]=5;
	      nrep21=1;
	    }
	  //printf("initial sumZ1, sum1, sum, n1, n0, n1+n0=% lf %lf %lf %d %d %d \n",sumZ1, sum1, sum, nZ0, nZ1, nZ0+nZ1, n0, n1, n0+n1);
	  if (*met==0)
	    {
	      lambda2[i*I+0]=sumZrep2[i]/nZrep21[i];
	      lambda2[i*I+1]=0;
	      mu2[i*I+1]=sumrep2[i]/nrep21;
	      mu2[i*I+0]=0;
	      phi2[i*I+1]=1;
	      phi2[i*I+0]=0;
	    }
	  if (*met==1)
	    {
	      lambda2[i*I+0]=sumZrep2[i]/nZrep21[i];//sum0/n0;
	      lambda2[i*I+1]=sumrep2[i]/nrep21; 
	      mu2[i*I+0]=0;
	      mu2[i*I+1]=0;
	      phi2[i*I+0]=0;
	      phi2[i*I+1]=0;
	    }
	  if (*met==2)
	    {
	      lambda2[i*I+0]=0;//sum0/n0;
	      lambda2[i*I+1]=0; 
	      mu2[i*I+0]=sumZrep2[i]/nZrep21[i];
	      mu2[i*I+1]=sumrep2[i]/nrep21;
	      phi2[i*I+0]=1;
	      phi2[i*I+1]=1;
	    }
	  pi2[i]=(nZrep21[i]+0.0)/(nZrep20[i]+0.0);
	}
      qrep2[0]=0.1; //prob(X_s=1|X_{s-1}=0)
      qrep2[1]=0.9; //prob(X_s=1|X_{s-1}=1)
    }

  // set initial values for parameter of non-replicates of first condition(protein)
  if(nsp1>0)
    {
      for (i=0;i<nsp1;i++)
	{   
	  i1=i+nsr1;
	  sumZnrep1[i]=0.0;
	  sumnrep1[i]=0.0;
	  nnrep11[i]=0;
	  nnrep10[i]=0;
	  nZnrep11[i]=0;
	  nZnrep10[i]=0;
	  for (s=0;s<S;s++)
	    { 
	      if(dataset1[i1*S+s]>cp)
		{ 
		  stateXnrep1[i*S+s]=1;
		  nnrep11[i]=nnrep11[i]++;
		  sumnrep1[i]=sumnrep1[i]+dataset1[i1*S+s];
		  stateZ1[i1*S+s]=2;
		  indexY1[i1*S+s]=0;
		}
	      else
		{
		  stateXnrep1[i*S+s]=0;
		  nnrep10[i]=nnrep10[i]++;
		  if (dataset1[i1*S+s]==0)
		    {
		      indexY1[i1*S+s]=1;
		      stateZ1[i1*S+s]=0;
		      nZnrep10[i]=nZnrep10[i]++;
		    }
		  else
		    {
		      indexY1[i1*S+s]=0;
		      stateZ1[i1*S+s]=1;
		      sumZnrep1[i]=sumZnrep1[i]+dataset1[i1*S+s];
		      nZnrep11[i]=nZnrep11[i]++;
		    }
		}
	    }
	  if (sumnrep1[i]==0|nnrep11[i]==0)
	    {
	      sumnrep1[i]=5;
	      nnrep11[i]=1;
	    }
	  //printf("initial sumZ1, sum1, sum, n1, n0, n1+n0=% lf %lf %lf %d %d %d \n",sumZ1, sum1, sum, nZ0, nZ1, nZ0+nZ1, n0, n1, n0+n1);
	  if (*met==0)
	    {
	      lambda1[i1*I+0]=sumZnrep1[i]/nZnrep11[i];
	      lambda1[i1*I+1]=0;
	      mu1[i1*I+1]=sumnrep1[i]/nnrep11[i];
	      mu1[i1*I+0]=0;
	      phi1[i1*I+1]=1;
	      phi1[i1*I+0]=0;
	    }
	  if (*met==1)
	    {
	      lambda1[i1*I+0]=sumZnrep1[i]/nZnrep11[i];//sum0/n0;
	      lambda1[i1*I+1]=sumnrep1[i]/nnrep11[i]; 
	      mu1[i1*I+0]=0;
	      mu1[i1*I+1]=0;
	      phi1[i1*I+0]=0;
	      phi1[i1*I+1]=0;
	      //printf("initial lambda %lf %lf\n", lambda[0],lambda[1]); 
	    }
	  if (*met==2)
	    {
	      lambda1[i1*I+0]=0;//sum0/n0;
	      lambda1[i1*I+1]=0; 
	      mu1[i1*I+0]=sumZnrep1[i]/nZnrep11[i];
	      mu1[i1*I+1]=sumnrep1[i]/nnrep11[i];
	      phi1[i1*I+0]=1;
	      phi1[i1*I+1]=1;
	      //printf("initial mu and pi %lf %lf %lf\n", mu[0],mu[1], phi[0], phi[1]);
	    }
	  pi1[i1]=(nZnrep11[i]+0.0)/(nZnrep10[i]+0.0);
	  //printf("initial nZ and pi %d %d %lf\n", nZ1, nZ0, pi);
	  qnrep1[i*I+0]=0.1; //prob(X_s=1|X_{s-1}=0)
	  qnrep1[i*I+1]=0.9; //prob(X_s=1|X_{s-1}=1)
	}
    }
  // set initial values for parameter of non-replicates of second condition(protein)
  if (nsp2>0)
    {
      for (i=0;i<nsp2;i++)
	{   
	  i2=i+nsr2;
	  sumZnrep2[i]=0.0;
	  sumnrep2[i]=0.0;
	  nnrep21[i]=0;
	  nnrep20[i]=0;
	  nZnrep21[i]=0;
	  nZnrep20[i]=0;
	  for (s=0;s<S;s++)
	    { 
	      if(dataset2[i2*S+s]>cp)
		{ 
		  stateXnrep2[i*S+s]=1;
		  nnrep21[i]=nnrep21[i]++;
		  sumnrep2[i]=sumnrep2[i]+dataset2[i2*S+s];
		  stateZ2[i2*S+s]=2;
		  indexY2[i2*S+s]=0;
		}
	      else
		{
		  stateXnrep2[i*S+s]=0;
		  nnrep20[i]=nnrep20[i]++;
		  if (dataset2[i2*S+s]==0)
		    {
		      indexY2[i2*S+s]=1;
		      stateZ2[i2*S+s]=0;
		      nZnrep20[i]=nZnrep20[i]++;
		    }
		  else
		    {
		      indexY2[i2*S+s]=0;
		      stateZ2[i2*S+s]=1;
		      sumZnrep2[i]=sumZnrep2[i]+dataset2[i2*S+s];
		      nZnrep21[i]=nZnrep21[i]++;
		    }
		}
	    }
	  if (sumnrep2[i]==0|nnrep21[i]==0)
	    {
	      sumnrep2[i]=5;
	      nnrep21[i]=1;
	    }
	  //printf("initial sumZ1, sum1, sum, n1, n0, n1+n0=% lf %lf %lf %d %d %d \n",sumZ1, sum1, sum, nZ0, nZ1, nZ0+nZ1, n0, n1, n0+n1);
	  if (*met==0)
	    {
	      lambda2[i2*I+0]=sumZnrep2[i]/nZnrep21[i];
	      lambda2[i2*I+1]=0;
	      mu2[i2*I+1]=sumnrep2[i]/nnrep21[i];
	      mu2[i2*I+0]=0;
	      phi2[i2*I+1]=1;
	      phi2[i2*I+0]=0;
	    }
	  if (*met==1)
	    {
	      lambda2[i2*I+0]=sumZnrep2[i]/nZnrep21[i];//sum0/n0;
	      lambda2[i2*I+1]=sumnrep2[i]/nnrep21[i]; 
	      mu2[i2*I+0]=0;
	      mu2[i2*I+1]=0;
	      phi2[i2*I+0]=0;
	      phi2[i2*I+1]=0;
	      //printf("initial lambda %lf %lf\n", lambda[0],lambda[1]); 
	    }
	  if (*met==2)
	    {
	      lambda2[i2*I+0]=0;//sum0/n0;
	      lambda2[i2*I+1]=0; 
	      mu2[i2*I+0]=sumZnrep2[i]/nZnrep21[i];
	      mu2[i2*I+1]=sumnrep2[i]/nnrep21[i];
	      phi2[i2*I+0]=1;
	      phi2[i2*I+1]=1;
	      //printf("initial mu and pi %lf %lf %lf\n", mu[0],mu[1], phi[0], phi[1]);
	    }
	  pi2[i2]=(nZnrep21[i]+0.0)/(nZnrep20[i]+0.0);
	  //printf("initial nZ and pi %d %d %lf\n", nZ1, nZ0, pi);
	  qnrep2[i*I+0]=0.1; //prob(X_s=1|X_{s-1}=0)
	  qnrep2[i*I+1]=0.9; //prob(X_s=1|X_{s-1}=1)
	}
    }

  // get variances for Metrolis-Hastling method
  double *sdmu10, *sdmu11, *sdphi10, *sdphi11,*rate1;
  sdmu10=(double*)malloc(sizeof(double*)*nr1);
  sdmu11=(double*)malloc(sizeof(double*)*nr1);
  sdphi10=(double*)malloc(sizeof(double*)*nr1); 
  sdphi11=(double*)malloc(sizeof(double*)*nr1); 
  int nr14=4*nr1;
  rate1=(double*)malloc(sizeof(double*)*nr14);
  for (i=0;i<nr1;i++)
    {
      sdmu11[i]=*(var1+i*4+0);
      sdphi11[i]=*(var1+i*4+1);
      sdmu10[i]=*(var1+i*4+2);
      sdphi10[i]=*(var1+i*4+3);
      rate1[i*4+0]=0.0;
      rate1[i*4+1]=0.0;
      rate1[i*4+2]=0.0;
      rate1[i*4+3]=0.0;
      //Rprintf("var of first condition= %lf %lf %lf %lf \n", sdmu11[i], sdphi11[i], sdmu10[i], sdphi10[i]);
    }
  double *sdmu20, *sdmu21, *sdphi20, *sdphi21,*rate2;
  sdmu20=(double*)malloc(sizeof(double*)*nr2);
  sdmu21=(double*)malloc(sizeof(double*)*nr2);
  sdphi20=(double*)malloc(sizeof(double*)*nr2); 
  sdphi21=(double*)malloc(sizeof(double*)*nr2); 
  int nr24=4*nr2;
  rate2=(double*)malloc(sizeof(double*)*nr24);
  for (i=0;i<nr2;i++)
    {
      sdmu21[i]=*(var2+i*4+0);
      sdphi21[i]=*(var2+i*4+1);
      sdmu20[i]=*(var2+i*4+2);
      sdphi20[i]=*(var2+i*4+3);
      rate2[i*4+0]=0.0;
      rate2[i*4+1]=0.0;
      rate2[i*4+2]=0.0;
      rate2[i*4+3]=0.0;
      //Rprintf("var of second condition= %lf %lf %lf %lf \n", sdmu21[i], sdphi21[i], sdmu20[i], sdphi20[i]);
    }	    

  double  sdq0rep1, *sdq0nrep1, sdq0rep2, *sdq0nrep2, *rate;
  sdq0nrep1=(double*)malloc(sizeof(double*)*nsp1);
  sdq0nrep2=(double*)malloc(sizeof(double*)*nsp2);
  int nrate=indexnsr1+indexnsr2+nsp1+nsp2+1;
  rate=(double*)malloc(sizeof(double*)*nrate);
  for (i=0;i<nrate;i++)
    {
      rate[i]=0;
    }
 // Rprintf("var of transition probability =\n");
  if (nsr1>0)
    {
      sdq0rep1=*(var+0);
     //Rprintf("%lf ", sdq0rep1);
    }
  if (nsp1>0)
    {
      for (i=0;i<nsp1;i++)
	{
	  sdq0nrep1[i]=*(var+indexnsr1+i);
	  //Rprintf("%lf ", sdq0nrep1[i]);
	}
    }
  if (nsr2>0)
    {
      sdq0rep2=*(var+indexnsr1+nsp1+0);
      //Rprintf("%lf ", sdq0rep2);
    }
  if (nsp2>0)
    {
      for (i=0;i<nsp2;i++)
	{
	  sdq0nrep2[i]=*(var+indexnsr1+nsp1+indexnsr2+i);
	  //Rprintf("%lf ", sdq0nrep2[i]);
	}
    }
  double sdc=*(var+nrate-1);
  //Rprintf("%lf\n", sdc);

  double propmu0, propmu1, propphi0, propphi1, Cmu0, Cmu1, logAmu0, logAmu1, Cphi0, Cphi1, logAphi0, logAphi1, sumgammaphi0, sumgammaphi1;
  double propc, Cc, logAc, logAc1, logAc2, propq0, Cq0, logAq0;


  j1=0;
  int z1=0;
  GetRNGstate();
  for(z=0;z<*N;z++)//show results for every 10th iterations
    {
      // 0. Sample the common protion c for transition probabilities. This c is common for two conditions and we sample c by Metropolis-Hasting method
      // 0.1 sample proposal value for c
      U1=runif(0,1);
      temp=U1*pnorm(c/sdc, 0, 1, 1, 0)+pnorm(-c/sdc, 0, 1, 1, 0);
      propc=qnorm(temp, 0, 1, 1, 0)*sdc+c;
      logAc1=0;
      logAc2=0;

      //**************************************************sample states and parameters for replicates of first condition********************************************************************************
      if (nsr1>0)
	{
	  // 1. sample X and Z
	  // 1.1 sample stateX and stateZ for replicates
	  Q[0][0]=1-qrep1[0];
	  Q[0][1]=qrep1[0];
	  Q[1][0]=1-qrep1[1];
	  Q[1][1]=qrep1[1];
	  n[0][0]=0;
	  n[0][1]=0;
	  n[1][0]=0;
	  n[1][1]=0;
	  nrep10=0;
	  nrep11=0;

	  for (s=0; s<S; s++)
	    {
	      for (j=0;j<nsr1;j++)
		{
		  if (*met==0)// poissonNB
		    {
		      logpY1[0][j*S+s]=dpois(dataset1[j*S+s], lambda1[j*I+0], 1);//use R function dpois(y, lambda, logistrue=1);
		      pY1[0][j*S+s]=dpois(dataset1[j*S+s], lambda1[j*I+0], 0);//use R function dpois(y, lambda, logisfalse=0);
		      logpY1[1][j*S+s]=dnbinom(dataset1[j*S+s], phi1[j*I+1], phi1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
		      pY1[1][j*S+s]=dnbinom(dataset1[j*S+s], phi1[j*I+1], phi1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
		    }
		  if (*met==1)// poisson
		    {
		      for (i=0;i<I;i++)
			{
			  logpY1[i][j*S+s]=dpois(dataset1[j*S+s], lambda1[j*I+i], 1);//use R function dpois(y, lambda, logistrue=1);
			  pY1[i][j*S+s]=dpois(dataset1[j*S+s], lambda1[j*I+i], 0);//use R function dpois(y, lambda, logisfalse=0);
			}
		    }
		  if (*met==2)// NB
		    {
		      for (i=0;i<I;i++)
			{
			  logpY1[i][j*S+s]=dnbinom(dataset1[j*S+s], phi1[j*I+i], phi1[j*I+i]/(mu1[j*I+i]+phi1[j*I+i]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
			  pY1[i][j*S+s]=dnbinom(dataset1[j*S+s], phi1[j*I+i], phi1[j*I+i]/(mu1[j*I+i]+phi1[j*I+i]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
			}
		    }
		}
	    }
	  pXcYrep1[0][0]=1.0;
	  pXcYrep1[1][0]=1.0;
	  for (j=0;j<nsr1;j++)
	    {
	      pXcYrep1[0][0]=pXcYrep1[0][0]*((1-pi1[j])*indexY1[j*S+0]+pi1[j]*indexY1[j*S+0]*pY1[0][j*S+0]+(1-indexY1[j*S+0])*pY1[0][j*S+0]);
	      pXcYrep1[1][0]=pXcYrep1[1][0]*pY1[1][j*S+0];
	    }
	  pXcYrep1[0][0]=pXcYrep1[0][0]*Q[0][stateXrep1[1]];
	  pXcYrep1[1][0]=pXcYrep1[1][0]*Q[1][stateXrep1[1]];
	  probXrep1[0]=1/(1+exp(log(pXcYrep1[0][0])-log(pXcYrep1[1][0])));//posterior distribution of P(X_1=1)
	  sumprobXrep1[0]=sumprobXrep1[0]+probXrep1[0];
	  U=runif(0,1);      
	  if(U<probXrep1[0])
	    {
	      stateXrep1[0]=1;
	      nrep11=nrep11++;//To get sum_s I(X_s=1)
	    }
	  else
	    {
	      stateXrep1[0]=0;
	      nrep10=nrep10++;
	    }
	  for (s=1; s<S-1; s++)
	    {
	      pXcYrep1[0][s]=1.0;
	      pXcYrep1[1][s]=1.0;
	      for (j=0;j<nsr1;j++)
		{
		  pXcYrep1[0][s]=pXcYrep1[0][s]*((1-pi1[j])*indexY1[j*S+s]+pi1[j]*indexY1[j*S+s]*pY1[0][j*S+s]+(1-indexY1[j*S+s])*pY1[0][j*S+s]);
		  pXcYrep1[1][s]=pXcYrep1[1][s]*pY1[1][j*S+s];
		}
	      pXcYrep1[0][s]=pXcYrep1[0][s]*Q[stateXrep1[s-1]][0]*Q[0][stateXrep1[s+1]];
	      pXcYrep1[1][s]=pXcYrep1[1][s]*Q[stateXrep1[s-1]][1]*Q[1][stateXrep1[s+1]];
	      probXrep1[s]=1/(1+exp(log(pXcYrep1[0][s])-log(pXcYrep1[1][s])));//posterior distribution of P(X_1=1)
	      sumprobXrep1[s]=sumprobXrep1[s]+probXrep1[s];
	      U=runif(0,1);
	      if(U<probXrep1[s])
		{
		  stateXrep1[s]=1;
		  nrep11=nrep11++;//To get sum_s I(X_s=1)
		}     
	      else
		{
		  stateXrep1[s]=0; 
		  nrep10=nrep10++;
		}
	     
	    }
	  pXcYrep1[0][S-1]=1.0;
	  pXcYrep1[1][S-1]=1.0;
	  for (j=0;j<nsr1;j++)
	    {
	      pXcYrep1[0][S-1]=pXcYrep1[0][S-1]*((1-pi1[j])*indexY1[j*S+S-1]+pi1[j]*indexY1[j*S+S-1]*pY1[0][j*S+S-1]+(1-indexY1[j*S+S-1])*pY1[0][j*S+S-1]);
	      pXcYrep1[1][S-1]=pXcYrep1[1][S-1]*pY1[1][j*S+S-1];
	    }
	  pXcYrep1[0][S-1]=pXcYrep1[0][S-1]*Q[stateXrep1[S-2]][0];
	  pXcYrep1[1][S-1]=pXcYrep1[1][S-1]*Q[stateXrep1[S-2]][1];
	  probXrep1[S-1]=1/(1+exp(log(pXcYrep1[0][S-1])-log(pXcYrep1[1][S-1])));//posterior distribution of P(X_1=1)
	  sumprobXrep1[S-1]=sumprobXrep1[S-1]+probXrep1[S-1];
	  U=runif(0,1);
	  if(U<probXrep1[S-1])
	    {
	      stateXrep1[S-1]=1;
	      nrep11=nrep11++;//To get sum_s I(X_s=1)
	    }     
	  else
	    {	  
	      stateXrep1[S-1]=0;
	      nrep10=nrep10++;
	    }
      
	  // 1.2 given stateX we sample stateZ
	  for (j=0;j<nsr1;j++)
	    {
	      sumZrep1[j]=0.0;
	      sumrep1[j]=0.0;
	      nZrep10[j]=0;
	      nZrep11[j]=0;
	      if (*met<2)
		{
		  probZ1[j]=1/(1+exp(log(1-pi1[j])-log(pi1[j])+lambda1[j*I+0]));
		}
	      if (*met==2)
		{
		  probZ1[j]=1/(1+exp(log(1-pi1[j])-log(pi1[j])-phi1[j*I+0]*log(phi1[j*I+0]/(phi1[j*I+0]+mu1[j*I+0]))));
		}
	      //Rprintf("probZ= %lf\n", probZ);
	      for (s=0;s<S;s++)
		{
		  if (stateXrep1[s]==0)
		    {
		      if (dataset1[j*S+s]>0)
			{
			  stateZ1[j*S+s]=1;
			  sumZrep1[j]=sumZrep1[j]+dataset1[j*S+s];
			  nZrep11[j]=nZrep11[j]++;
			}
		      else
			{
			  U=runif(0,1);
			  if (U<probZ1[j])
			    {
			      stateZ1[j*S+s]=1;
			      sumZrep1[j]=sumZrep1[j]+dataset1[j*S+s]; //To get sum_s Y_sI(X_s=0, Z_s=1)
			      nZrep11[j]=nZrep11[j]++;//To get sum_s I(X_s=0, Z_s=1)
			    }
			  else
			    {
			      stateZ1[j*S+s]=0;
			      nZrep10[j]=nZrep10[j]++;//To get sum_s I(X_s=0, Z_s=0)
			    }
			}
		    }
		  else
		    {		  
		      sumrep1[j]=sumrep1[j]+dataset1[j*S+s]; //To get sum_s Y_sI(X_s=1)
		      stateZ1[j*S+s]=2;
		    }
		}
	    }
	  
	  //*************************************************************************************************************************	
	  // iteration step 2: sample the parameters from posterior distribution.
	  // ***************************** sampling parameters for replicates of first condition
	  // 2.1.1 count the pair
	  for (s=0;s<S-1;s++)
	    {
	      sum2X[s]=stateXrep1[s]+stateXrep1[s+1]*2;
	      if(sum2X[s]==3)
		{
		  n[1][1]++;
		}
	      else if(sum2X[s]==2)
		{
		  n[0][1]++;
		}
	      else if(sum2X[s]==1)
		{
		  n[1][0]++;
		}
	      else
		{
		  n[0][0]++;
		}
	    }
	  //Rprintf("n=%d %d %d %d\n", n[1][1], n[1][0], n[0][1], n[0][0]);
	  
	  // 2.1.2 sample transition probability q0
	  U1=runif(0,1);
	  temp=U1*(pnorm((1-qrep1[0])/sdq0rep1, 0, 1, 1, 0)-pnorm(-qrep1[0]/sdq0rep1, 0, 1, 1, 0))+pnorm(-qrep1[0]/sdq0rep1, 0, 1, 1, 0);//U1*pnorm(qrep1[0]/sdq0rep1, 0, 1, 1, 0)+pnorm(-qrep1[0]/sdq0rep1, 0, 1, 1, 0);
	  propq0=qnorm(temp, 0, 1, 1, 0)*sdq0rep1+qrep1[0];	
	  // 2.1.2. Calculate the acceptance ratio for q0
	  Cq0=log(pnorm((1-qrep1[0])/sdq0rep1, 0, 1, 1, 0)-pnorm(-qrep1[0]/sdq0rep1,0,1,1,0))-log(pnorm((1-propq0)/sdq0rep1, 0, 1, 1, 0)-pnorm(-propq0/sdq0rep1,0,1,1,0));//pnorm(qrep1[0]/sdq0rep1,0,1,1,1)-pnorm(propq0/sdq0rep1,0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)//truncated Normal(mu, sigma, 0, 1)
	  //logAq0=n[1][1]*(log(c-propq0)-log(c-qrep1[0]))+(n[1][0]+n[0][1]+AQrep1-1)*(log(propq0)-log(qrep1[0]))+(n[0][0]+BQrep1-1)*(log(1-propq0)-log(1-qrep1[0]))+Cq0;
	  logAq0=n[1][1]*(log(1-propq0*c)-log(1-qrep1[0]*c))+(n[1][0]+n[0][1]+AQrep1-1)*(log(propq0)-log(qrep1[0]))+(n[0][0]+BQrep1-1)*(log(1-propq0)-log(1-qrep1[0]))+Cq0;// 2.1.3. sample q0
	  U=runif(0,1);
	  if(log(U)<logAq0)
	    {
	      qrep1[0]=propq0;
	      rate[0]=rate[0]+1.0;
	    }
	  // 2.1.4. claculate q1 by use the q1=1-q0/c
	  qrep1[1]=1-qrep1[0]*c;
	  // 2.1.4. claculate logAc1 and logAc2 which is used later when sample the common protion c.  
	  //logAc1=logAc1+stateXrep1[0]*log(1/(1+1/propc))+(1-stateXrep1[0])*log(1-1/(1+1/propc))+n[1][1]*log(1-qrep1[0]/propc)+n[1][0]*log(qrep1[0]/propc);
	  //logAc2=logAc2+stateXrep1[0]*log(1/(1+1/c))+(1-stateXrep1[0])*log(1-1/(1+1/c))+n[1][1]*log(1-qrep1[0]/c)+n[1][0]*log(qrep1[0]/c);
	  logAc1=logAc1+stateXrep1[0]*log(1/(1+propc))+(1-stateXrep1[0])*log(1-1/(1+propc))+n[1][1]*log(1-qrep1[0]*propc)+n[1][0]*log(qrep1[0]*propc);
	  logAc2=logAc2+stateXrep1[0]*log(1/(1+c))+(1-stateXrep1[0])*log(1-1/(1+c))+n[1][1]*log(1-qrep1[0]*c)+n[1][0]*log(qrep1[0]*c);
	  for (j=0;j<nsr1;j++)
	    {
	      // 2.2 sample inner mixture portion pi
	      a=Api1[j]+nZrep11[j];
	      b=Bpi1[j]+nZrep10[j];
	      pi1[j]=rbeta(a,b);
	      //printf("pi= %lf\n", pi[j]);
	      
	      // 2.3. sample lambda0 and mu1,phi1 when method="ZIP&NB"
	      if (*met==0)
		{
		  // 2.3.1. sample lambda0
		  a=Alambda1[j*I+0]+sumZrep1[j];
		  b=1/(Blambda1[j*I+0]+nZrep11[j]);
		  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda[0], sumZ1, Blambda[0], nZ1);
		  lambda1[j*I+0]=rgamma(a,b);
		  
		  // 2.3.2. sample mu1
		  U1=runif(0,1);
		  temp=U1*pnorm(mu1[j*I+1]/sdmu11[j], 0, 1, 1, 0)+pnorm(-mu1[j*I+1]/sdmu11[j], 0, 1, 1, 0);
		  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu11[j]+mu1[j*I+1];	  
		  Cmu1=pnorm(mu1[j*I+1]/sdmu11[j],0,1,1,1)-pnorm(propmu1/sdmu11[j],0,1,1,1);
		  logAmu1=sumrep1[j]*(log(propmu1/(propmu1+phi1[j*I+1]))-log(mu1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1])))+phi1[j*I+1]*nrep11*(log(mu1[j*I+1]+phi1[j*I+1])-log(propmu1+phi1[j*I+1]))+(Amu1[j*I+1]-1.0)*log(propmu1/mu1[j*I+1])-Bmu1[j*I+1]*(propmu1-mu1[j*I+1])+Cmu1;
		  U=runif(0,1);
		  if(log(U)<logAmu1)
		    {
		      mu1[j*I+1]=propmu1;
		      rate1[j*4+0]=rate1[j*4+0]+1.0;
		    }
		  // 2.3.3. sample phi1
		  U1=runif(0,1);
		  temp=U1*pnorm(phi1[j*I+1]/sdphi11[j], 0, 1, 1, 0)+pnorm(-phi1[j*I+1]/sdphi11[j], 0, 1, 1, 0);
		  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi11[j]+phi1[j*I+1];	
		  
		  Cphi1=pnorm(phi1[j*I+1]/sdphi11[j],0,1,1,1)-pnorm(propphi1/sdphi11[j],0,1,1,1);
		  sumgammaphi1=0;
		  for (s=0;s<S;s++)
		    {
		      if (stateXrep1[s]==1&dataset1[j*S+s]>0)
			{
			  for (j2=0;j2<dataset1[j*S+s];j2++)
			    {
			      sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi1[j*I+1]+j2);
			    }
			}
		    }
		  logAphi1=sumgammaphi1+sumrep1[j]*(log(mu1[j*I+1]+phi1[j*I+1])-log(mu1[j*I+1]+propphi1))+nrep11*(propphi1*log(propphi1/(mu1[j*I+1]+propphi1))-phi1[j*I+1]*log(phi1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1])))+(Aphi1[j*I+1]-1.0)*log(propphi1/phi1[j*I+1])-Bphi1[j*I+1]*(propphi1-phi1[j*I+1])+Cphi1;
		  U=runif(0,1);
		  if(log(U)<logAphi1)
		    {
		      phi1[j*I+1]=propphi1;
		      rate1[j*4+1]=rate1[j*4+1]+1.0;
		    }
		}
	  	  
	      // 2.4. sample lambda when method="poisson"(*met=1)
	      if (*met==1)
		{
		  a=Alambda1[j*I+0]+sumZrep1[j];
		  b=1/(Blambda1[j*I+0]+nZrep11[j]);
		  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda[0], sumZ1, Blambda[0], nZ1);
		  lambda1[j*I+0]=rgamma(a,b);
		  a=Alambda1[j*I+1]+sumrep1[j];
		  b=1/(Blambda1[j*I+1]+nrep11);
		  lambda1[j*I+1]=rgamma(a,b);
		  //Rprintf("sample lambda1 %lf %lf %lf %d\n", Alambda[1], sum1, Blambda[1], n1);
		  if(lambda1[j*I+0]>lambda1[j*I+1])
		    {
		      temp=lambda1[j*I+1];
		      lambda1[j*I+1]=lambda1[j*I+0];
		      lambda1[j*I+0]=temp;
		      //Rprintf("lambda= %lf %lf\n", lambda[0], lambda[1]);
		    }
		  //Rprintf("lambda= %lf %lf\n", lambda[0], lambda[1]);
		}
	  
	      // 2.5. sample mu and phi when method="NB" (*met=2)
	      if (*met==2) // sample hiarchique parameters phi_0 and phi_1 of Gamma distribution by Metropolis-Hastling method  
		{
		  // 2.5.1.Sample proposal mu from truncated normal distributions (postive tail)
		  U1=runif(0,1);
		  temp=U1*pnorm(mu1[j*I+0]/sdmu10[j], 0, 1, 1, 0)+pnorm(-mu1[j*I+0]/sdmu10[j], 0, 1, 1, 0);
		  propmu0=qnorm(temp, 0, 1, 1, 0)*sdmu10[j]+mu1[j*I+0];	
		  U1=runif(0,1);
		  temp=U1*pnorm(mu1[j*I+1]/sdmu11[j], 0, 1, 1, 0)+pnorm(-mu1[j*I+1]/sdmu11[j], 0, 1, 1, 0);
		  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu11[j]+mu1[j*I+1];	
		  
		  // 2.5.2. Calculate the acceptance ratio for mu
		  Cmu0=pnorm(mu1[j*I+0]/sdmu10[j],0,1,1,1)-pnorm(propmu0/sdmu10[j],0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
		  Cmu1=pnorm(mu1[j*I+1]/sdmu11[j],0,1,1,1)-pnorm(propmu1/sdmu11[j],0,1,1,1);
		  logAmu0=sumZrep1[j]*(log(propmu0/(propmu0+phi1[j*I+0]))-log(mu1[j*I+0]/(mu1[j*I+0]+phi1[j*I+0])))+phi1[j*I+0]*nZrep11[j]*(log(mu1[j*I+0]+phi1[j*I+0])-log(propmu0+phi1[j*I+0]))+(Amu1[j*I+0]-1.0)*log(propmu0/mu1[j*I+0])-Bmu1[j*I+0]*(propmu0-mu1[j*I+0])+Cmu0;
		  logAmu1=sumrep1[j]*(log(propmu1/(propmu1+phi1[j*I+1]))-log(mu1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1])))+phi1[j*I+1]*nrep11*(log(mu1[j*I+1]+phi1[j*I+1])-log(propmu1+phi1[j*I+1]))+(Amu1[j*I+1]-1.0)*log(propmu1/mu1[j*I+1])-Bmu1[j*I+1]*(propmu1-mu1[j*I+1])+Cmu1;
		  // 2.5.3. sample mu
		  U=runif(0,1);
		  if(log(U)<logAmu0)
		    {
		      mu1[j*I+0]=propmu0;
		      rate1[j*4+2]=rate1[j*4+2]+1.0;
		    }
		  U=runif(0,1);
		  if(log(U)<logAmu1)
		    {
		      mu1[j*I+1]=propmu1;
		      rate1[j*4+0]=rate1[j*4+0]+1.0;
		    }
		  
	      // 2.5.4. Sample proposal phi from truncated normal distribution (postive tail)
		  U1=runif(0,1);
		  temp=U1*pnorm(phi1[j*I+0]/sdphi10[j], 0, 1, 1, 0)+pnorm(-phi1[j*I+0]/sdphi10[j], 0, 1, 1, 0);
		  propphi0=qnorm(temp, 0, 1, 1, 0)*sdphi10[j]+phi1[j*I+0];	
		  U1=runif(0,1);
		  temp=U1*pnorm(phi1[j*I+1]/sdphi11[j], 0, 1, 1, 0)+pnorm(-phi1[j*I+1]/sdphi11[j], 0, 1, 1, 0);
		  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi11[j]+phi1[j*I+1];	
	      
		  // 2.5.5. Calculate the acceptance ratio for phi
		  Cphi0=pnorm(phi1[j*I+0]/sdphi10[j],0,1,1,1)-pnorm(propphi0/sdphi10[j],0,1,1,1);// log(PHI(phi^t/sd_phi))-log(PHI(phi'/sd_phi))##using R function pnorm(x, mu, phi, iflowtail=1, iflog=1)
		  Cphi1=pnorm(phi1[j*I+1]/sdphi11[j],0,1,1,1)-pnorm(propphi1/sdphi11[j],0,1,1,1);
		  sumgammaphi0=0;
		  sumgammaphi1=0;
		  for (s=0;s<S;s++)
		    {
		      if (stateXrep1[s]==0&stateZ1[j*S+s]==1&dataset1[j*S+s]>0)
			{
			  for (j2=0;j2<dataset1[j*S+s];j2++)
			    {
			      sumgammaphi0=sumgammaphi0+log(propphi0+j2)-log(phi1[j*I+0]+j2);
			    }
			}
		      if (stateXrep1[s]==1&dataset1[j*S+s]>0)
			{
			  for (j2=0;j2<dataset1[j*S+s];j2++)
			    {
			      sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi1[j*I+1]+j2);
			    }
			}	
		    }
		  logAphi0=sumgammaphi0+sumZrep1[j]*(log(mu1[j*I+0]+phi1[j*I+0])-log(mu1[j*I+0]+propphi0))+nZrep11[j]*(propphi0*log(propphi0/(mu1[j*I+0]+propphi0))-phi1[j*I+0]*log(phi1[j*I+0]/(mu1[j*I+0]+phi1[j*I+0])))+(Aphi1[j*I+0]-1.0)*log(propphi0/phi1[j*I+0])-Bphi1[j*I+0]*(propphi0-phi1[j*I+0])+Cphi0;
		  logAphi1=sumgammaphi1+sumrep1[j]*(log(mu1[j*I+1]+phi1[j*I+1])-log(mu1[j*I+1]+propphi1))+nrep11*(propphi1*log(propphi1/(mu1[j*I+1]+propphi1))-phi1[j*I+1]*log(phi1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1])))+(Aphi1[j*I+1]-1.0)*log(propphi1/phi1[j*I+1])-Bphi1[j*I+1]*(propphi1-phi1[j*I+1])+Cphi1;
		  // 2.5.6. sample phi
		  U=runif(0,1);
		  if(log(U)<logAphi0)
		    {
		      phi1[j*I+0]=propphi0;
		      rate1[j*4+3]=rate1[j*4+3]+1.0;
		    }
		  U=runif(0,1);
		  if(log(U)<logAphi1)
		    {
		      phi1[j*I+1]=propphi1;
		      rate1[j*4+1]=rate1[j*4+1]+1.0;
		    }
		  //Rprintf("mu,phi= %lf %lf %lf %lf\n",mu[j*I+1], phi[j*I+1], mu[j*I+0],  phi[j*I+0]);
		  if(mu1[j*I+0]>mu1[j*I+1])
		    {
		      //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu[0], mu[1], phi[0], phi[1]);
		      temp=phi1[j*I+1];
		      phi1[j*I+1]=phi1[j*I+0];
		      phi1[j*I+0]=temp;
		      temp=mu1[j*I+1];
		      mu1[j*I+1]=mu1[j*I+0];
		      mu1[j*I+0]=temp;
		    }
		  //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu[0], mu[1], phi[0], phi[1]);
		}
	    }
	}



      //************************************************************** sample states and parameters for replicates of second condition******************************************************************
      if (nsr2>1)
	{
	  Q[0][0]=1-qrep2[0];
	  Q[0][1]=qrep2[0];
	  Q[1][0]=1-qrep2[1];
	  Q[1][1]=qrep2[1];
	  n[0][0]=0;
	  n[0][1]=0;
	  n[1][0]=0;
	  n[1][1]=0;
	  nrep20=0;
	  nrep21=0;

	  for (s=0; s<S; s++)
	    {
	      for (j=0;j<nsr2;j++)
		{
		  if (*met==0)// poissonNB
		    {
		      logpY2[0][j*S+s]=dpois(dataset2[j*S+s], lambda2[j*I+0], 1);//use R function dpois(y, lambda, logistrue=1);
		      pY2[0][j*S+s]=dpois(dataset2[j*S+s], lambda2[j*I+0], 0);//use R function dpois(y, lambda, logisfalse=0);
		      logpY2[1][j*S+s]=dnbinom(dataset2[j*S+s], phi2[j*I+1], phi2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
		      pY2[1][j*S+s]=dnbinom(dataset2[j*S+s], phi2[j*I+1], phi2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
		    }
		  if (*met==1)// poisson
		    {
		      for (i=0;i<I;i++)
			{
			  logpY2[i][j*S+s]=dpois(dataset2[j*S+s], lambda2[j*I+i], 1);//use R function dpois(y, lambda, logistrue=1);
			  pY2[i][j*S+s]=dpois(dataset2[j*S+s], lambda2[j*I+i], 0);//use R function dpois(y, lambda, logisfalse=0);
			}
		    }
		  if (*met==2)// NB
		    {
		      for (i=0;i<I;i++)
			{
			  logpY2[i][j*S+s]=dnbinom(dataset2[j*S+s], phi2[j*I+i], phi2[j*I+i]/(mu2[j*I+i]+phi2[j*I+i]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
			  pY2[i][j*S+s]=dnbinom(dataset2[j*S+s], phi2[j*I+i], phi2[j*I+i]/(mu2[j*I+i]+phi2[j*I+i]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
			}
		    }
		}
	    }
	  pXcYrep2[0][0]=1.0;
	  pXcYrep2[1][0]=1.0;
	  for (j=0;j<nsr2;j++)
	    {
	      pXcYrep2[0][0]=pXcYrep2[0][0]*((1-pi2[j])*indexY2[j*S+0]+pi2[j]*indexY2[j*S+0]*pY2[0][j*S+0]+(1-indexY2[j*S+0])*pY2[0][j*S+0]);
	      pXcYrep2[1][0]=pXcYrep2[1][0]*pY2[1][j*S+0];
	    }
	  pXcYrep2[0][0]=pXcYrep2[0][0]*Q[0][stateXrep2[1]];
	  pXcYrep2[1][0]=pXcYrep2[1][0]*Q[1][stateXrep2[1]];
	  probXrep2[0]=1/(1+exp(log(pXcYrep2[0][0])-log(pXcYrep2[1][0])));//posterior distribution of P(X_1=1)
	  sumprobXrep2[0]=sumprobXrep2[0]+probXrep2[0];
	  U=runif(0,1);      
	  if(U<probXrep2[0])
	    {
	      stateXrep2[0]=1;
	      nrep21=nrep21++;//To get sum_s I(X_s=1)
	    }
	  else
	    {
	      stateXrep2[0]=0;
	      nrep20=nrep20++;
	    }
	  for (s=1; s<S-1; s++)
	    {
	      pXcYrep2[0][s]=1.0;
	      pXcYrep2[1][s]=1.0;
	      for (j=0;j<nsr2;j++)
		{
		  pXcYrep2[0][s]=pXcYrep2[0][s]*((1-pi2[j])*indexY2[j*S+s]+pi2[j]*indexY2[j*S+s]*pY2[0][j*S+s]+(1-indexY2[j*S+s])*pY2[0][j*S+s]);
		  pXcYrep2[1][s]=pXcYrep2[1][s]*pY2[1][j*S+s];
		}
	      pXcYrep2[0][s]=pXcYrep2[0][s]*Q[stateXrep2[s-1]][0]*Q[0][stateXrep2[s+1]];
	      pXcYrep2[1][s]=pXcYrep2[1][s]*Q[stateXrep2[s-1]][1]*Q[1][stateXrep2[s+1]];
	      probXrep2[s]=1/(1+exp(log(pXcYrep2[0][s])-log(pXcYrep2[1][s])));//posterior distribution of P(X_1=1)
	      sumprobXrep2[s]=sumprobXrep2[s]+probXrep2[s];
	      U=runif(0,1);
	      if(U<probXrep2[s])
		{
		  stateXrep2[s]=1;
		  nrep21=nrep21++;//To get sum_s I(X_s=1)
		}     
	      else
		{
		  stateXrep2[s]=0; 
		  nrep20=nrep20++;
		}
	      
	    }
	  pXcYrep2[0][S-1]=1.0;
	  pXcYrep2[1][S-1]=1.0;
	  for (j=0;j<nsr2;j++)
	    {
	      pXcYrep2[0][S-1]=pXcYrep2[0][S-1]*((1-pi2[j])*indexY2[j*S+S-1]+pi2[j]*indexY2[j*S+S-1]*pY2[0][j*S+S-1]+(1-indexY2[j*S+S-1])*pY2[0][j*S+S-1]);
	      pXcYrep2[1][S-1]=pXcYrep2[1][S-1]*pY2[1][j*S+S-1];
	    }
	  pXcYrep2[0][S-1]=pXcYrep2[0][S-1]*Q[stateXrep2[S-2]][0];
	  pXcYrep2[1][S-1]=pXcYrep2[1][S-1]*Q[stateXrep2[S-2]][1];
	  probXrep2[S-1]=1/(1+exp(log(pXcYrep2[0][S-1])-log(pXcYrep2[1][S-1])));//posterior distribution of P(X_1=1)
	  sumprobXrep2[S-1]=sumprobXrep2[S-1]+probXrep2[S-1];
	  U=runif(0,1);
	  if(U<probXrep2[S-1])
	    {
	      stateXrep2[S-1]=1;
	      nrep21=nrep21++;//To get sum_s I(X_s=1)
	    }     
	  else
	    {	  
	      stateXrep2[S-1]=0;
	      nrep20=nrep20++;
	    }
	  
	  // 1.2 given stateX we sample stateZ
	  for (j=0;j<nsr2;j++)
	    {
	      sumZrep2[j]=0.0;
	      sumrep2[j]=0.0;
	      nZrep20[j]=0;
	      nZrep21[j]=0;
	      if (*met<2)
		{
		  probZ2[j]=1/(1+exp(log(1-pi2[j])-log(pi2[j])+lambda2[j*I+0]));
		}
	      if (*met==2)
		{
		  probZ2[j]=1/(1+exp(log(1-pi2[j])-log(pi2[j])-phi2[j*I+0]*log(phi2[j*I+0]/(phi2[j*I+0]+mu2[j*I+0]))));
		}
	      //Rprintf("probZ= %lf\n", probZ);
	      for (s=0;s<S;s++)
		{
		  if (stateXrep2[s]==0)
		    {
		      if (dataset2[j*S+s]>0)
			{
			  stateZ2[j*S+s]=1;
			  sumZrep2[j]=sumZrep2[j]+dataset2[j*S+s];
			  nZrep21[j]=nZrep21[j]++;
			}
		      else
			{
			  U=runif(0,1);
			  if (U<probZ2[j])
			    {
			      stateZ2[j*S+s]=1;
			      sumZrep2[j]=sumZrep2[j]+dataset2[j*S+s]; //To get sum_s Y_sI(X_s=0, Z_s=1)
			      nZrep21[j]=nZrep21[j]++;//To get sum_s I(X_s=0, Z_s=1)
			    }
			  else
			    {
			      stateZ2[j*S+s]=0;
			      nZrep20[j]=nZrep20[j]++;//To get sum_s I(X_s=0, Z_s=0)
			    }
			}
		    }
		  else
		    {		  
		      sumrep2[j]=sumrep2[j]+dataset2[j*S+s]; //To get sum_s Y_sI(X_s=1)
		      stateZ2[j*S+s]=2;
		    }
		}
	    }
	  
	  //*************************************************************************************************************************	
          // iteration step 2: sample the parameters from posterior distribution.
	  // ***************************** sampling parameters for replicates of second condition
	  // 2.1.1 count the pair
	  for (s=0;s<S-1;s++)
	    {
	      sum2X[s]=stateXrep2[s]+stateXrep2[s+1]*2;
	      if(sum2X[s]==3)
		{
		  n[1][1]++;
		}
	      else if(sum2X[s]==2)
		{
		  n[0][1]++;
		}
	      else if(sum2X[s]==1)
		{
		  n[1][0]++;
		}
	      else
		{
		  n[0][0]++;
		}
	    }
	  //Rprintf("n=%d %d %d %d\n", n[1][1], n[1][0], n[0][1], n[0][0]);
	  
	  // 2.1.2 sample transition probability q0
	  U1=runif(0,1);
	  temp=U1*(pnorm((1-qrep2[0])/sdq0rep2, 0, 1, 1, 0)-pnorm(-qrep2[0]/sdq0rep2, 0, 1, 1, 0))+pnorm(-qrep2[0]/sdq0rep2, 0, 1, 1, 0);//U1*pnorm(qrep2[0]/sdq0rep2, 0, 1, 1, 0)+pnorm(-qrep2[0]/sdq0rep2, 0, 1, 1, 0);
	  propq0=qnorm(temp, 0, 1, 1, 0)*sdq0rep2+qrep2[0];	
	  // 2.1.2. Calculate the acceptance ratio for q0
	  Cq0=log(pnorm((1-qrep2[0])/sdq0rep2, 0, 1, 1, 0)-pnorm(-qrep2[0]/sdq0rep2,0,1,1,0))-log(pnorm((1-propq0)/sdq0rep2, 0, 1, 1, 0)-pnorm(-propq0/sdq0rep2,0,1,1,0));//pnorm(qrep2[0]/sdq0rep2,0,1,1,1)-pnorm(propq0/sdq0rep2,0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
	  logAq0=n[1][1]*(log(1-propq0*c)-log(1-qrep2[0]*c))+(n[1][0]+n[0][1]+AQrep2-1)*(log(propq0)-log(qrep2[0]))+(n[0][0]+BQrep2-1)*(log(1-propq0)-log(1-qrep2[0]))+Cq0;
	  // 2.1.3. sample q0
	  U=runif(0,1);
	  if(log(U)<logAq0)
	    {
	      qrep2[0]=propq0;
	      rate[indexnsr1+nsp1+0]=rate[indexnsr1+nsp1+0]+1.0;
	    }
	  // 2.1.4. claculate q1 by use the q1=1-q0/c
	  qrep2[1]=1-qrep2[0]*c;
	  // 2.1.4. claculate logAc1 and logAc2 which is used later when sample the common protion c.  
	  //logAc1=logAc1+stateXrep2[0]*log(1/(1+1/propc))+(1-stateXrep2[0])*log(1-1/(1+1/propc))+n[1][1]*log(1-qrep2[0]/propc)+n[1][0]*log(qrep2[0]/propc);
	  //logAc2=logAc2+stateXrep2[0]*log(1/(1+1/c))+(1-stateXrep2[0])*log(1-1/(1+1/c))+n[1][1]*log(1-qrep2[0]/c)+n[1][0]*log(qrep2[0]/c);
	  logAc1=logAc1+stateXrep2[0]*log(1/(1+propc))+(1-stateXrep2[0])*log(1-1/(1+propc))+n[1][1]*log(1-qrep2[0]*propc)+n[1][0]*log(qrep2[0]*propc);
	  logAc2=logAc2+stateXrep2[0]*log(1/(1+c))+(1-stateXrep2[0])*log(1-1/(1+c))+n[1][1]*log(1-qrep2[0]*c)+n[1][0]*log(qrep2[0]*c);
	  
	  for (j=0;j<nsr2;j++)
	    {
	      // 2.2 sample inner mixture portion pi
	      a=Api2[j]+nZrep21[j];
	      b=Bpi2[j]+nZrep20[j];
	      pi2[j]=rbeta(a,b);
	      //printf("pi= %lf\n", pi2[j]);
	      
	      // 2.3. sample lambda0 and mu1,phi1 when method="ZIP&NB"
	      if (*met==0)
		{
		  // 2.3.1. sample lambda0
		  a=Alambda2[j*I+0]+sumZrep2[j];
		  b=1/(Blambda2[j*I+0]+nZrep21[j]);
		  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda[0], sumZrep2, Blambda2[0], nZrep2);
		  lambda2[j*I+0]=rgamma(a,b);
		  
		  // 2.3.2. sample mu1
		  U1=runif(0,1);
		  temp=U1*pnorm(mu2[j*I+1]/sdmu21[j], 0, 1, 1, 0)+pnorm(-mu2[j*I+1]/sdmu21[j], 0, 1, 1, 0);
		  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu21[j]+mu2[j*I+1];	  
		  Cmu1=pnorm(mu2[j*I+1]/sdmu21[j],0,1,1,1)-pnorm(propmu1/sdmu21[j],0,1,1,1);
		  logAmu1=sumrep2[j]*(log(propmu1/(propmu1+phi2[j*I+1]))-log(mu2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1])))+phi2[j*I+1]*nrep21*(log(mu2[j*I+1]+phi2[j*I+1])-log(propmu1+phi2[j*I+1]))+(Amu2[j*I+1]-1.0)*log(propmu1/mu2[j*I+1])-Bmu2[j*I+1]*(propmu1-mu2[j*I+1])+Cmu1;
		  U=runif(0,1);
		  if(log(U)<logAmu1)
		    {
		      mu2[j*I+1]=propmu1;
		      rate2[j*4+0]=rate2[j*4+0]+1.0;
		    }
		  // 2.3.3. sample phi1
		  U1=runif(0,1);
		  temp=U1*pnorm(phi2[j*I+1]/sdphi21[j], 0, 1, 1, 0)+pnorm(-phi2[j*I+1]/sdphi21[j], 0, 1, 1, 0);
		  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi21[j]+phi2[j*I+1];	
		  
		  Cphi1=pnorm(phi2[j*I+1]/sdphi21[j],0,1,1,1)-pnorm(propphi1/sdphi21[j],0,1,1,1);
		  sumgammaphi1=0;
		  for (s=0;s<S;s++)
		    {
		      if (stateXrep2[s]==1&dataset2[j*S+s]>0)
			{
			  for (j2=0;j2<dataset2[j*S+s];j2++)
			    {
			      sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi2[j*I+1]+j2);
			    }
			}
		    }
		  logAphi1=sumgammaphi1+sumrep2[j]*(log(mu2[j*I+1]+phi2[j*I+1])-log(mu2[j*I+1]+propphi1))+nrep21*(propphi1*log(propphi1/(mu2[j*I+1]+propphi1))-phi2[j*I+1]*log(phi2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1])))+(Aphi2[j*I+1]-1.0)*log(propphi1/phi2[j*I+1])-Bphi2[j*I+1]*(propphi1-phi2[j*I+1])+Cphi1;
		  U=runif(0,1);
		  if(log(U)<logAphi1)
		    {
		      phi2[j*I+1]=propphi1;
		      rate2[j*4+1]=rate2[j*4+1]+1.0;
		    }
		}

	      // 2.4. sample lambda when method="poisson"(*met=1)
	      if (*met==1)
		{
		  a=Alambda2[j*I+0]+sumZrep2[j];
		  b=1/(Blambda2[j*I+0]+nZrep21[j]);
		  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda[0], sumZ1, Blambda[0], nZ1);
		  lambda2[j*I+0]=rgamma(a,b);
		  a=Alambda2[j*I+1]+sumrep2[j];
		  b=1/(Blambda2[j*I+1]+nrep21);
		  lambda2[j*I+1]=rgamma(a,b);
		  //Rprintf("sample lambda1 %lf %lf %lf %d\n", Alambda[1], sum1, Blambda[1], n1);
		  if(lambda2[j*I+0]>lambda2[j*I+1])
		    {
		      temp=lambda2[j*I+1];
		      lambda2[j*I+1]=lambda2[j*I+0];
		      lambda2[j*I+0]=temp;
		      //Rprintf("lambda= %lf %lf\n", lambda2[0], lambda2[1]);
		    }
		  //Rprintf("lambda= %lf %lf\n", lambda2[0], lambda2[1]);
		}
	  
	      // 2.5. sample mu and phi when method="NB" (*met=2)
	      if (*met==2) // sample hiarchique parameters phi_0 and phi_1 of Gamma distribution by Metropolis-Hastling method  
		{
		  // 2.5.1.Sample proposal mu from truncated normal distributions (postive tail)
		  U1=runif(0,1);
		  temp=U1*pnorm(mu2[j*I+0]/sdmu20[j], 0, 1, 1, 0)+pnorm(-mu2[j*I+0]/sdmu20[j], 0, 1, 1, 0);
		  propmu0=qnorm(temp, 0, 1, 1, 0)*sdmu20[j]+mu2[j*I+0];	
		  U1=runif(0,1);
		  temp=U1*pnorm(mu2[j*I+1]/sdmu21[j], 0, 1, 1, 0)+pnorm(-mu2[j*I+1]/sdmu21[j], 0, 1, 1, 0);
		  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu21[j]+mu2[j*I+1];	
		  
		  // 2.5.2. Calculate the acceptance ratio for mu
		  Cmu0=pnorm(mu2[j*I+0]/sdmu20[j],0,1,1,1)-pnorm(propmu0/sdmu20[j],0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
		  Cmu1=pnorm(mu2[j*I+1]/sdmu21[j],0,1,1,1)-pnorm(propmu1/sdmu21[j],0,1,1,1);
		  logAmu0=sumZrep2[j]*(log(propmu0/(propmu0+phi2[j*I+0]))-log(mu2[j*I+0]/(mu2[j*I+0]+phi2[j*I+0])))+phi2[j*I+0]*nZrep21[j]*(log(mu2[j*I+0]+phi2[j*I+0])-log(propmu0+phi2[j*I+0]))+(Amu2[j*I+0]-1.0)*log(propmu0/mu2[j*I+0])-Bmu2[j*I+0]*(propmu0-mu2[j*I+0])+Cmu0;
		  logAmu1=sumrep2[j]*(log(propmu1/(propmu1+phi2[j*I+1]))-log(mu2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1])))+phi2[j*I+1]*nrep21*(log(mu2[j*I+1]+phi2[j*I+1])-log(propmu1+phi2[j*I+1]))+(Amu2[j*I+1]-1.0)*log(propmu1/mu2[j*I+1])-Bmu2[j*I+1]*(propmu1-mu2[j*I+1])+Cmu1;
		  // 2.5.3. sample mu
		  U=runif(0,1);
		  if(log(U)<logAmu0)
		    {
		      mu2[j*I+0]=propmu0;
		      rate2[j*4+2]=rate2[j*4+2]+1.0;
		    }
		  U=runif(0,1);
		  if(log(U)<logAmu1)
		    {
		      mu2[j*I+1]=propmu1;
		      rate2[j*4+0]=rate2[j*4+0]+1.0;
		    }
	      
		  // 2.5.4. Sample proposal phi from truncated normal distribution (postive tail)
		  U1=runif(0,1);
		  temp=U1*pnorm(phi2[j*I+0]/sdphi20[j], 0, 1, 1, 0)+pnorm(-phi2[j*I+0]/sdphi20[j], 0, 1, 1, 0);
		  propphi0=qnorm(temp, 0, 1, 1, 0)*sdphi20[j]+phi2[j*I+0];	
		  U1=runif(0,1);
		  temp=U1*pnorm(phi2[j*I+1]/sdphi21[j], 0, 1, 1, 0)+pnorm(-phi2[j*I+1]/sdphi21[j], 0, 1, 1, 0);
		  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi21[j]+phi2[j*I+1];	
		  
		  // 2.5.5. Calculate the acceptance ratio for phi
		  Cphi0=pnorm(phi2[j*I+0]/sdphi20[j],0,1,1,1)-pnorm(propphi0/sdphi20[j],0,1,1,1);// log(PHI(phi^t/sd_phi))-log(PHI(phi'/sd_phi))##using R function pnorm(x, mu, phi, iflowtail=1, iflog=1)
		  Cphi1=pnorm(phi2[j*I+1]/sdphi21[j],0,1,1,1)-pnorm(propphi1/sdphi21[j],0,1,1,1);
		  sumgammaphi0=0;
		  sumgammaphi1=0;
		  for (s=0;s<S;s++)
		    {
		      if (stateXrep2[s]==0&stateZ2[j*S+s]==1&dataset2[j*S+s]>0)
			{
			  for (j2=0;j2<dataset2[j*S+s];j2++)
			    {
			      sumgammaphi0=sumgammaphi0+log(propphi0+j2)-log(phi2[j*I+0]+j2);
			    }
			}
		      if (stateXrep2[s]==1&dataset2[j*S+s]>0)
			{
			  for (j2=0;j2<dataset2[j*S+s];j2++)
			    {
			      sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi2[j*I+1]+j2);
			    }
			}	
		    }
		  logAphi0=sumgammaphi0+sumZrep2[j]*(log(mu2[j*I+0]+phi2[j*I+0])-log(mu2[j*I+0]+propphi0))+nZrep21[j]*(propphi0*log(propphi0/(mu2[j*I+0]+propphi0))-phi2[j*I+0]*log(phi2[j*I+0]/(mu2[j*I+0]+phi2[j*I+0])))+(Aphi2[j*I+0]-1.0)*log(propphi0/phi2[j*I+0])-Bphi2[j*I+0]*(propphi0-phi2[j*I+0])+Cphi0;
		  logAphi1=sumgammaphi1+sumrep2[j]*(log(mu2[j*I+1]+phi2[j*I+1])-log(mu2[j*I+1]+propphi1))+nrep21*(propphi1*log(propphi1/(mu2[j*I+1]+propphi1))-phi2[j*I+1]*log(phi2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1])))+(Aphi2[j*I+1]-1.0)*log(propphi1/phi2[j*I+1])-Bphi2[j*I+1]*(propphi1-phi2[j*I+1])+Cphi1;
		  // 2.5.6. sample phi
		  U=runif(0,1);
		  if(log(U)<logAphi0)
		    {
		      phi2[j*I+0]=propphi0;
		      rate2[j*4+3]=rate2[j*4+3]+1.0;
		    }
		  U=runif(0,1);
		  if(log(U)<logAphi1)
		    {
		      phi2[j*I+1]=propphi1;
		      rate2[j*4+1]=rate2[j*4+1]+1.0;
		    }
		  //Rprintf("mu,phi= %lf %lf %lf %lf\n",mu2[j*I+1], phi2[j*I+1], mu2[j*I+0],  phi2[j*I+0]);
		  if(mu2[j*I+0]>mu2[j*I+1])
		    {
		      //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu[0], mu[1], phi[0], phi[1]);
		      temp=phi2[j*I+1];
		      phi2[j*I+1]=phi2[j*I+0];
		      phi2[j*I+0]=temp;
		      temp=mu2[j*I+1];
		      mu2[j*I+1]=mu2[j*I+0];
		      mu2[j*I+0]=temp;
		    }
		  //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu2[0], mu2[1], phi2[0], phi2[1]);
		}
	    }
	}


      /////////////////////////////////////////////start check

      //********************************************************* sample states and parameters for nonreplicates of first condition**************************************************************************
      if (nsp1>0)
	{
	  for (j=nsr1;j<nr1;j++)
	    {	    
	      i1=j-nsr1;
	      n[0][0]=0;
	      n[0][1]=0;
	      n[1][0]=0;
	      n[1][1]=0;
	      Q[0][0]=1-qnrep1[i1*I+0];
	      Q[0][1]=qnrep1[i1*I+0];
	      Q[1][0]=1-qnrep1[i1*I+1];
	      Q[1][1]=qnrep1[i1*I+1];	 
	      nnrep10[i1]=0;
	      nnrep11[i1]=0;
	      sumZnrep1[i1]=0.0;
	      sumnrep1[i1]=0.0;
	      nZnrep10[i1]=0;
	      nZnrep11[i1]=0;
	      
	      for (s=0; s<S; s++)
		{
		  if (*met==0)// poissonNB
		    {
		      logpY1[0][j*S+s]=dpois(dataset1[j*S+s], lambda1[j*I+0], 1);//use R function dpois(y, lambda, logistrue=1);
		      pY1[0][j*S+s]=dpois(dataset1[j*S+s], lambda1[j*I+0], 0);//use R function dpois(y, lambda, logisfalse=0);
		      logpY1[1][j*S+s]=dnbinom(dataset1[j*S+s], phi1[j*I+1], phi1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
		      pY1[1][j*S+s]=dnbinom(dataset1[j*S+s], phi1[j*I+1], phi1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
		    }
		  if (*met==1)// poisson
		    {
		      for (i=0;i<I;i++)
			{
			  logpY1[i][j*S+s]=dpois(dataset1[j*S+s], lambda1[j*I+i], 1);//use R function dpois(y, lambda, logistrue=1);
			  pY1[i][j*S+s]=dpois(dataset1[j*S+s], lambda1[j*I+i], 0);//use R function dpois(y, lambda, logisfalse=0);
			}
		    }
		  if (*met==2)// NB
		    {
		      for (i=0;i<I;i++)
			{
			  logpY1[i][j*S+s]=dnbinom(dataset1[j*S+s], phi1[j*I+i], phi1[j*I+i]/(mu1[j*I+i]+phi1[j*I+i]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
			  pY1[i][j*S+s]=dnbinom(dataset1[j*S+s], phi1[j*I+i], phi1[j*I+i]/(mu1[j*I+i]+phi1[j*I+i]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
			}
		    }
		}
	      pXcYnrep1[0][i1*S+0]=((1-pi1[j])*indexY1[j*S+0]+pi1[j]*indexY1[j*S+0]*pY1[0][j*S+0]+(1-indexY1[j*S+0])*pY1[0][j*S+0])*Q[0][stateXnrep1[i1*S+1]];
	      pXcYnrep1[1][i1*S+0]=pY1[1][j*S+0]*Q[1][stateXnrep1[i1*S+1]];
	      probXnrep1[i1*S+0]=1/(1+exp(log(pXcYnrep1[0][i1*S+0])-log(pXcYnrep1[1][i1*S+0])));//posterior distribution of P(X_1=1)
	      sumprobXnrep1[i1*S+0]=sumprobXnrep1[i1*S+0]+probXnrep1[i1*S+0];
	      U=runif(0,1);      
	      if(U<probXnrep1[i1*S+0])
		{
		  stateXnrep1[i1*S+0]=1;
		  nnrep11[i1]=nnrep11[i1]++;//To get sum_s I(X_s=1)
		}
	      else
		{
		  stateXnrep1[i1*S+0]=0;
		  nnrep10[i1]=nnrep10[i1]++;
		}
	      for (s=1; s<S-1; s++)
		{
		  pXcYnrep1[0][i1*S+s]=((1-pi1[j])*indexY1[j*S+s]+pi1[j]*indexY1[j*S+s]*pY1[0][j*S+s]+(1-indexY1[j*S+s])*pY1[0][j*S+s])*Q[stateXnrep1[i1*S+s-1]][0]*Q[0][stateXnrep1[i1*S+s+1]];
		  pXcYnrep1[1][i1*S+s]=pY1[1][j*S+s]*Q[stateXnrep1[i1*S+s-1]][1]*Q[1][stateXnrep1[i1*S+s+1]];
		  probXnrep1[i1*S+s]=1/(1+exp(log(pXcYnrep1[0][i1*S+s])-log(pXcYnrep1[1][i1*S+s])));//posterior distribution of P(X_1=1)
		  sumprobXnrep1[i1*S+s]=sumprobXnrep1[i1*S+s]+probXnrep1[i1*S+s];
		  U=runif(0,1);
		  if(U<probXnrep1[i1*S+s])
		    {
		      stateXnrep1[i1*S+s]=1;
		      nnrep11[i1]=nnrep11[i1]++;//To get sum_s I(X_s=1)
		    }     
		  else
		    {
		      stateXnrep1[i1*S+s]=0; 
		      nnrep10[i1]=nnrep10[i1]++;
		    }
		}
	      pXcYnrep1[0][i1*S+S-1]=((1-pi1[j])*indexY1[j*S+S-1]+pi1[j]*indexY1[j*S+S-1]*pY1[0][j*S+S-1]+(1-indexY1[j*S+S-1])*pY1[0][j*S+S-1])*Q[stateXnrep1[i1*S+S-2]][0];
	      pXcYnrep1[1][i1*S+S-1]=pY1[1][j*S+S-1]*Q[stateXnrep1[i1*S+S-2]][1];
	      probXnrep1[i1*S+S-1]=1/(1+exp(log(pXcYnrep1[0][i1*S+S-1])-log(pXcYnrep1[1][i1*S+S-1])));//posterior distribution of P(X_1=1)
	      sumprobXnrep1[i1*S+S-1]=sumprobXnrep1[i1*S+S-1]+probXnrep1[i1*S+S-1];
	      U=runif(0,1);
	      if(U<probXnrep1[i1*S+S-1])
		{
		  stateXnrep1[i1*S+S-1]=1;
		  nnrep11[i1]=nnrep11[i1]++;//To get sum_s I(X_s=1)
		}     
	      else
		{	  
		  stateXnrep1[i1*S+S-1]=0;
		  nnrep10[i1]=nnrep10[i1]++;
		}
      
	      // 1.2 given stateX we sample stateZ
	      if (*met<2)
		{
		  probZ1[j]=1/(1+exp(log(1-pi1[j])-log(pi1[j])+lambda1[j*I+0]));
		}
	      if (*met==2)
		{
		  probZ1[j]=1/(1+exp(log(1-pi1[j])-log(pi1[j])-phi1[j*I+0]*log(phi1[j*I+0]/(phi1[j*I+0]+mu1[j*I+0]))));
		}
	      //Rprintf("probZ= %lf\n", probZ);
	      for (s=0;s<S;s++)
		{
		  if (stateXnrep1[i1*S+s]==0)
		    {
		      if (dataset1[i1*S+s]>0)
			{
			  stateZ1[j*S+s]=1;
			  sumZnrep1[i1]=sumZnrep1[i1]+dataset1[j*S+s];
			  nZnrep11[i1]=nZnrep11[i1]++;
			}
		      else
			{
			  U=runif(0,1);
			  if (U<probZ1[j])
			    {
			      stateZ1[j*S+s]=1;
			      sumZnrep1[i1]=sumZnrep1[i1]+dataset1[j*S+s]; //To get sum_s Y_sI(X_s=0, Z_s=1)
			      nZnrep11[i1]=nZnrep11[i1]++;//To get sum_s I(X_s=0, Z_s=1)
			    }
			  else
			    {
			      stateZ1[j*S+s]=0;
			      nZnrep10[i1]=nZnrep10[i1]++;//To get sum_s I(X_s=0, Z_s=0)
			    }
			}
		    }
		  else
		    {		  
		      sumnrep1[i1]=sumnrep1[i1]+dataset1[j*S+s]; //To get sum_s Y_sI(X_s=1)
		      stateZ1[j*S+s]=2;
		    }
		}
	      // *************************** sampling parameters for non-replicates of first condition
	      // 2.1.1 count the pair
	      for (s=0;s<S-1;s++)
		{
		  sum2X[s]=stateXnrep1[i1*S+s]+stateXnrep1[i1*S+s+1]*2;
		  if(sum2X[s]==3)
		    {
		      n[1][1]++;
		    }
		  else if(sum2X[s]==2)
		    {
		      n[0][1]++;
		    }
		  else if(sum2X[s]==1)
		    {
		      n[1][0]++;
		    }
		  else
		    {
		      n[0][0]++;
		    }
		}
	      //Rprintf("n=%d %d %d %d\n", n[1][1], n[1][0], n[0][1], n[0][0]);
	      
	      // 2.1.2 sample transition probability q_0, and later sample the common protion c. We sample q_0 by Metropolis-Hasting method
	      U1=runif(0,1);
	      temp==U1*(pnorm((1-qnrep1[i1*I+0])/sdq0nrep1[i1], 0, 1, 1, 0)-pnorm(-qnrep1[i1*I+0]/sdq0nrep1[i1], 0, 1, 1, 0))+pnorm(-qnrep1[i1*I+0]/sdq0nrep1[i1], 0, 1, 1, 0);//U1*pnorm(qnrep1[i1*I+0]/sdq0nrep1[i1], 0, 1, 1, 0)+pnorm(-qnrep1[i1*I+0]/sdq0nrep1[i1], 0, 1, 1, 0);
	      propq0=qnorm(temp, 0, 1, 1, 0)*sdq0nrep1[i1]+qnrep1[i1*I+0];		       
	      // 2.1.2. Calculate the acceptance ratio for q0
	      Cq0=log(pnorm((1-qnrep1[i1*I+0])/sdq0nrep1[i1], 0, 1, 1, 0)-pnorm(-qnrep1[i1*I+0]/sdq0nrep1[i1],0,1,1,0))-log(pnorm((1-propq0)/sdq0nrep1[i1], 0, 1, 1, 0)-pnorm(-propq0/sdq0nrep1[i1],0,1,1,0));//pnorm(qnrep1[i1*I+0]/sdq0nrep1[i1],0,1,1,1)-pnorm(propq0/sdq0nrep1[i1],0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)      
	      logAq0=n[1][1]*(log(1-propq0*c)-log(1-qnrep1[i1*I+0]*c))+(n[1][0]+n[0][1]+AQnrep1[i1]-1)*(log(propq0)-log(qnrep1[i1*I+0]))+(n[0][0]+BQnrep1[i1]-1)*(log(1-propq0)-log(1-qnrep1[i1*I+0]))+Cq0;// 2.1.3. sample q0
	      U=runif(0,1);
	      if(log(U)<logAq0)
		{
		  qnrep1[i1*I+0]=propq0;
		  rate[i1+indexnsr1]=rate[i1+indexnsr1]+1.0;
		}
	      // 2.1.4. claculate q1 by use the q1=1-q0/c
	      qnrep1[i1*I+1]=1-qnrep1[i1*I+0]*c;
	      // 2.1.4. claculate logAc1 and logAc2 which is used later when sample the common protion c.  
	      //logAc1=logAc1+stateXnrep1[i1*S+0]*log(1/(1+1/propc))+(1-stateXnrep1[i1*S+0])*log(1-1/(1+1/propc))+n[1][1]*log(1-qnrep1[i1*I+0]/propc)+n[1][0]*log(qnrep1[i1*I+0]/propc);
	      //logAc2=logAc2+stateXnrep1[i1*S+0]*log(1/(1+1/c))+(1-stateXnrep1[i1*S+0])*log(1-1/(1+1/c))+n[1][1]*log(1-qnrep1[i1*I+0]/c)+n[1][0]*log(qnrep1[i1*I+0]/c);
	      logAc1=logAc1+stateXnrep1[i1*S+0]*log(1/(1+propc))+(1-stateXnrep1[i1*S+0])*log(1-1/(1+propc))+n[1][1]*log(1-qnrep1[i1*I+0]*propc)+n[1][0]*log(qnrep1[i1*I+0]*propc);
	      logAc2=logAc2+stateXnrep1[i1*S+0]*log(1/(1+c))+(1-stateXnrep1[i1*S+0])*log(1-1/(1+c))+n[1][1]*log(1-qnrep1[i1*I+0]*c)+n[1][0]*log(qnrep1[i1*I+0]*c);
	      
	      // 2.2 sample inner mixture portion pi
	      a=Api1[j]+nZnrep11[i1];
	      b=Bpi1[j]+nZnrep10[i1];
	      pi1[j]=rbeta(a,b);
	      //printf("pi= %lf\n", pi[j]);
	      
	      // 2.3. sample lambda0 and mu1,phi1 when method="ZIP&NB"
	      if (*met==0)
		{
		  // 2.3.1. sample lambda0
		  a=Alambda1[j*I+0]+sumZnrep1[i1];
		  b=1/(Blambda1[j*I+0]+nZnrep11[i1]);
		  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda[0], sumZ1, Blambda[0], nZ1);
		  lambda1[j*I+0]=rgamma(a,b);
		  
		  // 2.3.2. sample mu1
		  U1=runif(0,1);
		  temp=U1*pnorm(mu1[j*I+1]/sdmu11[j], 0, 1, 1, 0)+pnorm(-mu1[j*I+1]/sdmu11[j], 0, 1, 1, 0);
		  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu11[j]+mu1[j*I+1];	  
		  Cmu1=pnorm(mu1[j*I+1]/sdmu11[j],0,1,1,1)-pnorm(propmu1/sdmu11[j],0,1,1,1);
		  logAmu1=sumnrep1[i1]*(log(propmu1/(propmu1+phi1[j*I+1]))-log(mu1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1])))+phi1[j*I+1]*nnrep11[i1]*(log(mu1[j*I+1]+phi1[j*I+1])-log(propmu1+phi1[j*I+1]))+(Amu1[j*I+1]-1.0)*log(propmu1/mu1[j*I+1])-Bmu1[j*I+1]*(propmu1-mu1[j*I+1])+Cmu1;
		  U=runif(0,1);
		  if(log(U)<logAmu1)
		    {
		      mu1[j*I+1]=propmu1;
		      rate1[j*4+0]=rate1[j*4+0]+1.0;
		    }
		  // 2.3.3. sample phi1
		  U1=runif(0,1);
		  temp=U1*pnorm(phi1[j*I+1]/sdphi11[j], 0, 1, 1, 0)+pnorm(-phi1[j*I+1]/sdphi11[j], 0, 1, 1, 0);
		  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi11[j]+phi1[j*I+1];	
		  
		  Cphi1=pnorm(phi1[j*I+1]/sdphi11[j],0,1,1,1)-pnorm(propphi1/sdphi11[j],0,1,1,1);
		  sumgammaphi1=0;
		  for (s=0;s<S;s++)
		    {
		      if (stateXnrep1[i1*S+s]==1&dataset1[j*S+s]>0)
			{
			  for (j2=0;j2<dataset1[j*S+s];j2++)
			    {
			      sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi1[j*I+1]+j2);
			    }
			}
		    }
		  logAphi1=sumgammaphi1+sumnrep1[i1]*(log(mu1[j*I+1]+phi1[j*I+1])-log(mu1[j*I+1]+propphi1))+nnrep11[i1]*(propphi1*log(propphi1/(mu1[j*I+1]+propphi1))-phi1[j*I+1]*log(phi1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1])))+(Aphi1[j*I+1]-1.0)*log(propphi1/phi1[j*I+1])-Bphi1[j*I+1]*(propphi1-phi1[j*I+1])+Cphi1;
		  U=runif(0,1);
		  if(log(U)<logAphi1)
		    {
		      phi1[j*I+1]=propphi1;
		      rate1[j*4+1]=rate1[j*4+1]+1.0;
		    }
		}

	      // 2.4. sample lambda when method="poisson"(*met=1)
	      if (*met==1)
		{
		  a=Alambda1[j*I+0]+sumZnrep1[i1];
		  b=1/(Blambda1[j*I+0]+nZnrep11[i1]);
		  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda1[0], sumZnrep1, Blambda1[0], nZnrep11);
		  lambda1[j*I+0]=rgamma(a,b);
		  a=Alambda1[j*I+1]+sumnrep1[i1];
		  b=1/(Blambda1[j*I+1]+nnrep11[i1]);
		  lambda1[j*I+1]=rgamma(a,b);
		  //Rprintf("sample lambda1 %lf %lf %lf %d\n", Alambda1[1], sumnrep1, Blambda1[1], nnrep11);
		  if(lambda1[j*I+0]>lambda1[j*I+1])
		    {
		      temp=lambda1[j*I+1];
		      lambda1[j*I+1]=lambda1[j*I+0];
		      lambda1[j*I+0]=temp;
		      //Rprintf("lambda= %lf %lf\n", lambda1[0], lambda1[1]);
		    }
		  //Rprintf("lambda= %lf %lf\n", lambda1[0], lambda1[1]);
		}
	      
	      // 2.5. sample mu and phi when method="NB" (*met=2)
	      if (*met==2) // sample hiarchique parameters phi_0 and phi_1 of Gamma distribution by Metropolis-Hastling method  
		{
		  // 2.5.1.Sample proposal mu from truncated normal distributions (postive tail)
		  U1=runif(0,1);
		  temp=U1*pnorm(mu1[j*I+0]/sdmu10[j], 0, 1, 1, 0)+pnorm(-mu1[j*I+0]/sdmu10[j], 0, 1, 1, 0);
		  propmu0=qnorm(temp, 0, 1, 1, 0)*sdmu10[j]+mu1[j*I+0];	
		  U1=runif(0,1);
		  temp=U1*pnorm(mu1[j*I+1]/sdmu11[j], 0, 1, 1, 0)+pnorm(-mu1[j*I+1]/sdmu11[j], 0, 1, 1, 0);
		  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu11[j]+mu1[j*I+1];	
		  
		  // 2.5.2. Calculate the acceptance ratio for mu
		  Cmu0=pnorm(mu1[j*I+0]/sdmu10[j],0,1,1,1)-pnorm(propmu0/sdmu10[j],0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
		  Cmu1=pnorm(mu1[j*I+1]/sdmu11[j],0,1,1,1)-pnorm(propmu1/sdmu11[j],0,1,1,1);
		  logAmu0=sumZnrep1[i1]*(log(propmu0/(propmu0+phi1[j*I+0]))-log(mu1[j*I+0]/(mu1[j*I+0]+phi1[j*I+0])))+phi1[j*I+0]*nZnrep11[i1]*(log(mu1[j*I+0]+phi1[j*I+0])-log(propmu0+phi1[j*I+0]))+(Amu1[j*I+0]-1.0)*log(propmu0/mu1[j*I+0])-Bmu1[j*I+0]*(propmu0-mu1[j*I+0])+Cmu0;
		  logAmu1=sumnrep1[i1]*(log(propmu1/(propmu1+phi1[j*I+1]))-log(mu1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1])))+phi1[j*I+1]*nnrep11[i1]*(log(mu1[j*I+1]+phi1[j*I+1])-log(propmu1+phi1[j*I+1]))+(Amu1[j*I+1]-1.0)*log(propmu1/mu1[j*I+1])-Bmu1[j*I+1]*(propmu1-mu1[j*I+1])+Cmu1;
		  // 2.5.3. sample mu
		  U=runif(0,1);
		  if(log(U)<logAmu0)
		    {
		      mu1[j*I+0]=propmu0;
		      rate1[j*4+2]=rate1[j*4+2]+1.0;
		    }
		  U=runif(0,1);
		  if(log(U)<logAmu1)
		    {
		      mu1[j*I+1]=propmu1;
		      rate1[j*4+0]=rate1[j*4+0]+1.0;
		    }
		  
		  // 2.5.4. Sample proposal phi from truncated normal distribution (postive tail)
		  U1=runif(0,1);
		  temp=U1*pnorm(phi1[j*I+0]/sdphi10[j], 0, 1, 1, 0)+pnorm(-phi1[j*I+0]/sdphi10[j], 0, 1, 1, 0);
		  propphi0=qnorm(temp, 0, 1, 1, 0)*sdphi10[j]+phi1[j*I+0];	
		  U1=runif(0,1);
		  temp=U1*pnorm(phi1[j*I+1]/sdphi11[j], 0, 1, 1, 0)+pnorm(-phi1[j*I+1]/sdphi11[j], 0, 1, 1, 0);
		  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi11[j]+phi1[j*I+1];	
		  
		  // 2.5.5. Calculate the acceptance ratio for phi
		  Cphi0=pnorm(phi1[j*I+0]/sdphi10[j],0,1,1,1)-pnorm(propphi0/sdphi10[j],0,1,1,1);// log(PHI(phi^t/sd_phi))-log(PHI(phi'/sd_phi))##using R function pnorm(x, mu, phi, iflowtail=1, iflog=1)
		  Cphi1=pnorm(phi1[j*I+1]/sdphi11[j],0,1,1,1)-pnorm(propphi1/sdphi11[j],0,1,1,1);
		  sumgammaphi0=0;
		  sumgammaphi1=0;
		  for (s=0;s<S;s++)
		    {
		      if (stateXnrep1[i1*S+s]==0&stateZ1[j*S+s]==1&dataset1[j*S+s]>0)
			{
			  for (j2=0;j2<dataset1[j*S+s];j2++)
			    {
			      sumgammaphi0=sumgammaphi0+log(propphi0+j2)-log(phi1[j*I+0]+j2);
			    }
			}
		      if (stateXnrep1[i1*S+s]==1&dataset1[j*S+s]>0)
			{
			  for (j2=0;j2<dataset1[j*S+s];j2++)
			    {
			      sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi1[j*I+1]+j2);
			    }
			}	
		    }
		  logAphi0=sumgammaphi0+sumZnrep1[i1]*(log(mu1[j*I+0]+phi1[j*I+0])-log(mu1[j*I+0]+propphi0))+nZnrep11[i1]*(propphi0*log(propphi0/(mu1[j*I+0]+propphi0))-phi1[j*I+0]*log(phi1[j*I+0]/(mu1[j*I+0]+phi1[j*I+0])))+(Aphi1[j*I+0]-1.0)*log(propphi0/phi1[j*I+0])-Bphi1[j*I+0]*(propphi0-phi1[j*I+0])+Cphi0;
		  logAphi1=sumgammaphi1+sumnrep1[i1]*(log(mu1[j*I+1]+phi1[j*I+1])-log(mu1[j*I+1]+propphi1))+nnrep11[i1]*(propphi1*log(propphi1/(mu1[j*I+1]+propphi1))-phi1[j*I+1]*log(phi1[j*I+1]/(mu1[j*I+1]+phi1[j*I+1])))+(Aphi1[j*I+1]-1.0)*log(propphi1/phi1[j*I+1])-Bphi1[j*I+1]*(propphi1-phi1[j*I+1])+Cphi1;
		  // 2.5.6. sample phi
		  U=runif(0,1);
		  if(log(U)<logAphi0)
		    {
		      phi1[j*I+0]=propphi0;
		      rate1[j*4+3]=rate1[j*4+3]+1.0;
		    }
		  U=runif(0,1);
		  if(log(U)<logAphi1)
		    {
		      phi1[j*I+1]=propphi1;
		      rate1[j*4+1]=rate1[j*4+1]+1.0;
		    }
		  //Rprintf("mu,phi= %lf %lf %lf %lf\n",mu1[j*I+1], phi1[j*I+1], mu1[j*I+0],  phi1[j*I+0]);
		  if(mu1[j*I+0]>mu1[j*I+1])
		    {
		      //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu1[0], mu1[1], phi1[0], phi1[1]);
		      temp=phi1[j*I+1];
		      phi1[j*I+1]=phi1[j*I+0];
		      phi1[j*I+0]=temp;
		      temp=mu1[j*I+1];
		      mu1[j*I+1]=mu1[j*I+0];
		      mu1[j*I+0]=temp;
		    }
		  //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu1[0], mu1[1], phi1[0], phi1[1]);
		}
	    }
	}

      //************************************************* sample nonreplicates for second condition*************************************************************************************************
      if (nsp2>0)
	{
	  for (j=nsr2;j<nr2;j++)
	    {	    
	      i2=j-nsr2;
	      n[0][0]=0;
	      n[0][1]=0;
	      n[1][0]=0;
	      n[1][1]=0;
	      Q[0][0]=1-qnrep2[i2*I+0];
	      Q[0][1]=qnrep2[i2*I+0];
	      Q[1][0]=1-qnrep2[i2*I+1];
	      Q[1][1]=qnrep2[i2*I+1];	 
	      nnrep10[i2]=0;
	      nnrep11[i2]=0;
	      sumZnrep1[i2]=0.0;
	      sumnrep1[i2]=0.0;
	      nZnrep10[i2]=0;
	      nZnrep11[i2]=0;

	      for (s=0; s<S; s++)
		{
		  if (*met==0)// poissonNB
		    {
		      logpY2[0][j*S+s]=dpois(dataset2[j*S+s], lambda2[j*I+0], 1);//use R function dpois(y, lambda, logistrue=1);
		      pY2[0][j*S+s]=dpois(dataset2[j*S+s], lambda2[j*I+0], 0);//use R function dpois(y, lambda, logisfalse=0);
		      logpY2[1][j*S+s]=dnbinom(dataset2[j*S+s], phi2[j*I+1], phi2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
		      pY2[1][j*S+s]=dnbinom(dataset2[j*S+s], phi2[j*I+1], phi2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
		    }
		  if (*met==1)// poisson
		    {
		      for (i=0;i<I;i++)
			{
			  logpY2[i][j*S+s]=dpois(dataset2[j*S+s], lambda2[j*I+i], 1);//use R function dpois(y, lambda, logistrue=1);
			  pY2[i][j*S+s]=dpois(dataset2[j*S+s], lambda2[j*I+i], 0);//use R function dpois(y, lambda, logisfalse=0);
			}
		    }
		  if (*met==2)// NB
		    {
		      for (i=0;i<I;i++)
			{
			  logpY2[i][j*S+s]=dnbinom(dataset2[j*S+s], phi2[j*I+i], phi2[j*I+i]/(mu2[j*I+i]+phi2[j*I+i]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
			  pY2[i][j*S+s]=dnbinom(dataset2[j*S+s], phi2[j*I+i], phi2[j*I+i]/(mu2[j*I+i]+phi2[j*I+i]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
			}
		    }
		}
	      pXcYnrep2[0][i2*S+0]=((1-pi2[j])*indexY2[j*S+0]+pi2[j]*indexY2[j*S+0]*pY2[0][j*S+0]+(1-indexY2[j*S+0])*pY2[0][j*S+0])*Q[0][stateXnrep2[i2*S+1]];
	      pXcYnrep2[1][i2*S+0]=pY2[1][j*S+0]*Q[1][stateXnrep2[i2*S+1]];
	      probXnrep2[i2*S+0]=1/(1+exp(log(pXcYnrep2[0][i2*S+0])-log(pXcYnrep2[1][i2*S+0])));//posterior distribution of P(X_1=1)
	      sumprobXnrep2[i2*S+0]=sumprobXnrep2[i2*S+0]+probXnrep2[i2*S+0];
	      U=runif(0,1);      
	      if(U<probXnrep2[i2*S+0])
		{
		  stateXnrep2[i2*S+0]=1;
		  nnrep21[i2]=nnrep21[i2]++;//To get sum_s I(X_s=1)
		}
	      else
		{
		  stateXnrep2[i2*S+0]=0;
		  nnrep20[i2]=nnrep20[i2]++;
		}
	      for (s=1; s<S-1; s++)
		{
		  pXcYnrep2[0][i2*S+s]=((1-pi2[j])*indexY2[j*S+s]+pi2[j]*indexY2[j*S+s]*pY2[0][j*S+s]+(1-indexY2[j*S+s])*pY2[0][j*S+s])*Q[stateXnrep2[i2*S+s-1]][0]*Q[0][stateXnrep2[i2*S+s+1]];
		  pXcYnrep2[1][i2*S+s]=pY2[1][j*S+s]*Q[stateXnrep2[i2*S+s-1]][1]*Q[1][stateXnrep2[i2*S+s+1]];
		  probXnrep2[i2*S+s]=1/(1+exp(log(pXcYnrep2[0][i2*S+s])-log(pXcYnrep2[1][i2*S+s])));//posterior distribution of P(X_1=1)
		  sumprobXnrep2[i2*S+s]=sumprobXnrep2[i2*S+s]+probXnrep2[i2*S+s];
		  U=runif(0,1);
		  if(U<probXnrep2[i2*S+s])
		    {
		      stateXnrep2[i2*S+s]=1;
		      nnrep21[i2]=nnrep21[i2]++;//To get sum_s I(X_s=1)
		    }     
		  else
		    {
		      stateXnrep2[i2*S+s]=0; 
		      nnrep20[i2]=nnrep20[i2]++;
		    }
		}
	      pXcYnrep2[0][i2*S+S-1]=((1-pi2[j])*indexY2[j*S+S-1]+pi2[j]*indexY2[j*S+S-1]*pY2[0][j*S+S-1]+(1-indexY2[j*S+S-1])*pY2[0][j*S+S-1])*Q[stateXnrep2[i2*S+S-2]][0];
	      pXcYnrep2[1][i2*S+S-1]=pY2[1][j*S+S-1]*Q[stateXnrep2[i2*S+S-2]][1];
	      probXnrep2[i2*S+S-1]=1/(1+exp(log(pXcYnrep2[0][i2*S+S-1])-log(pXcYnrep2[1][i2*S+S-1])));//posterior distribution of P(X_1=1)
	      sumprobXnrep2[i2*S+S-1]=sumprobXnrep2[i2*S+S-1]+probXnrep2[i2*S+S-1];
	      U=runif(0,1);
	      if(U<probXnrep2[i2*S+S-1])
		{
		  stateXnrep2[i2*S+S-1]=1;
		  nnrep21[i2]=nnrep21[i2]++;//To get sum_s I(X_s=1)
		}     
	      else
		{	  
		  stateXnrep2[i2*S+S-1]=0;
		  nnrep20[i2]=nnrep20[i2]++;
		}
	      
	      // 1.2 given stateX we sample stateZ
	      if (*met<2)
		{
		  probZ2[j]=1/(1+exp(log(1-pi2[j])-log(pi2[j])+lambda2[j*I+0]));
		}
	      if (*met==2)
		{
		  probZ2[j]=1/(1+exp(log(1-pi2[j])-log(pi2[j])-phi2[j*I+0]*log(phi2[j*I+0]/(phi2[j*I+0]+mu2[j*I+0]))));
		}
	      //Rprintf("probZ= %lf\n", probZ);
	      for (s=0;s<S;s++)
		{
		  if (stateXnrep2[i2*S+s]==0)
		    {
		      if (dataset2[i2*S+s]>0)
			{
			  stateZ2[j*S+s]=1;
			  sumZnrep2[i2]=sumZnrep2[i2]+dataset2[j*S+s];
			  nZnrep21[i2]=nZnrep21[i2]++;
			}
		      else
			{
			  U=runif(0,1);
			  if (U<probZ2[j])
			    {
			      stateZ2[j*S+s]=1;
			      sumZnrep2[i2]=sumZnrep2[i2]+dataset2[j*S+s]; //To get sum_s Y_sI(X_s=0, Z_s=1)
			      nZnrep21[i2]=nZnrep21[i2]++;//To get sum_s I(X_s=0, Z_s=1)
			    }
			  else
			    {
			      stateZ2[j*S+s]=0;
			      nZnrep20[i2]=nZnrep20[i2]++;//To get sum_s I(X_s=0, Z_s=0)
			    }
			}
		    }
		  else
		    {		  
		      sumnrep2[i2]=sumnrep2[i2]+dataset2[j*S+s]; //To get sum_s Y_sI(X_s=1)
		      stateZ2[j*S+s]=2;
		    }
		}

	      // *************************** sampling parameters for non-replicates of second condtion
	      // 2.1.1 count the pair
	      for (s=0;s<S-1;s++)
		{
		  sum2X[s]=stateXnrep2[i2*S+s]+stateXnrep2[i2*S+s+1]*2;
		  if(sum2X[s]==3)
		    {
		      n[1][1]++;
		    }
		  else if(sum2X[s]==2)
		    {
		      n[0][1]++;
		    }
		  else if(sum2X[s]==1)
		    {
		      n[1][0]++;
		    }
		  else
		    {
		      n[0][0]++;
		    }
		}
	      //Rprintf("n=%d %d %d %d\n", n[1][1], n[1][0], n[0][1], n[0][0]);
	      
	      // 2.1.2 sample transition probability q_0, and later sample the common protion c. We sample q_0 by Metropolis-Hasting method
	      U1=runif(0,1);
	      temp=U1*(pnorm((1-qnrep2[i2*I+0])/sdq0nrep2[i2], 0, 1, 1, 0)-pnorm(-qnrep2[i2*I+0]/sdq0nrep2[i2], 0, 1, 1, 0))+pnorm(-qnrep2[i2*I+0]/sdq0nrep2[i2], 0, 1, 1, 0);//U1*pnorm(qnrep2[i2*I+0]/sdq0nrep2[i2], 0, 1, 1, 0)+pnorm(-qnrep2[i2*I+0]/sdq0nrep2[i2], 0, 1, 1, 0);
	      propq0=qnorm(temp, 0, 1, 1, 0)*sdq0nrep2[i2]+qnrep2[i2*I+0];		       
	      // 2.1.2. Calculate the acceptance ratio for q0
	      Cq0=log(pnorm((1-qnrep2[i2*I+0])/sdq0nrep2[i2], 0, 1, 1, 0)-pnorm(-qnrep2[i2*I+0]/sdq0nrep2[i2],0,1,1,0))-log(pnorm((1-propq0)/sdq0nrep2[i2], 0, 1, 1, 0)-pnorm(-propq0/sdq0nrep2[i2],0,1,1,0));//pnorm(qnrep2[i2*I+0]/sdq0nrep2[i2],0,1,1,1)-pnorm(propq0/sdq0nrep2[i2],0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)#truncated normal
	      logAq0=n[1][1]*(log(1-propq0*c)-log(1-qnrep2[i2*I+0]*c))+(n[1][0]+n[0][1]+AQnrep2[i2]-1)*(log(propq0)-log(qnrep2[i2*I+0]))+(n[0][0]+BQnrep2[i2]-1)*(log(1-propq0)-log(1-qnrep2[i2*I+0]))+Cq0;
	      // 2.1.3. sample q0
	      U=runif(0,1);
	      if(log(U)<logAq0)
		{
		  qnrep2[i2*I+0]=propq0;
		  rate[indexnsr1+nsp1+indexnsr2+i2]=rate[indexnsr1+nsp1+indexnsr2+i2]+1.0;
		}
	      // 2.1.4. claculate q1 by use the q1=1-q0/c
	      qnrep2[i2*I+1]=1-qnrep2[i2*I+0]*c;
	      // 2.1.4. claculate logAc1 and logAc2 which is used later when sample the common protion c.  
	      //logAc1=logAc1+stateXnrep2[i2*S+0]*log(1/(1+1/propc))+(1-stateXnrep2[i2*S+0])*log(1-1/(1+1/propc))+n[1][1]*log(1-qnrep2[i2*I+0]/propc)+n[1][0]*log(qnrep2[i2*I+0]/propc);
	      //logAc2=logAc2+stateXnrep2[i2*S+0]*log(1/(1+1/c))+(1-stateXnrep2[i2*S+0])*log(1-1/(1+1/c))+n[1][1]*log(1-qnrep2[i2*I+0]/c)+n[1][0]*log(qnrep2[i2*I+0]/c);
	      logAc1=logAc1+stateXnrep2[i2*S+0]*log(1/(1+propc))+(1-stateXnrep2[i2*S+0])*log(1-1/(1+propc))+n[1][1]*log(1-qnrep2[i2*I+0]*propc)+n[1][0]*log(qnrep2[i2*I+0]*propc);
	      logAc2=logAc2+stateXnrep2[i2*S+0]*log(1/(1+c))+(1-stateXnrep2[i2*S+0])*log(1-1/(1+c))+n[1][1]*log(1-qnrep2[i2*I+0]*c)+n[1][0]*log(qnrep2[i2*I+0]*c);
	      
	      // 2.2 sample inner mixture portion pi
	      a=Api2[j]+nZnrep21[i2];
	      b=Bpi2[j]+nZnrep20[i2];
	      pi2[j]=rbeta(a,b);
	      //printf("pi= %lf\n", pi2[j]);

	      // 2.3. sample lambda0 and mu1,phi1 when method="ZIP&NB"
	      if (*met==0)
		{
		  // 2.3.1. sample lambda0
		  a=Alambda2[j*I+0]+sumZnrep2[i2];
		  b=1/(Blambda2[j*I+0]+nZnrep21[i2]);
		  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda2[0], sumZnrep21, Blambda2[0], nZnrep21);
		  lambda2[j*I+0]=rgamma(a,b);
		  
		  // 2.3.2. sample mu1
		  U1=runif(0,1);
		  temp=U1*pnorm(mu2[j*I+1]/sdmu21[j], 0, 1, 1, 0)+pnorm(-mu2[j*I+1]/sdmu21[j], 0, 1, 1, 0);
		  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu21[j]+mu2[j*I+1];	  
		  Cmu1=pnorm(mu2[j*I+1]/sdmu21[j],0,1,1,1)-pnorm(propmu1/sdmu21[j],0,1,1,1);
		  logAmu1=sumnrep2[i2]*(log(propmu1/(propmu1+phi2[j*I+1]))-log(mu2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1])))+phi2[j*I+1]*nnrep21[i2]*(log(mu2[j*I+1]+phi2[j*I+1])-log(propmu1+phi2[j*I+1]))+(Amu2[j*I+1]-1.0)*log(propmu1/mu2[j*I+1])-Bmu2[j*I+1]*(propmu1-mu2[j*I+1])+Cmu1;
		  U=runif(0,1);
		  if(log(U)<logAmu1)
		    {
		      mu2[j*I+1]=propmu1;
		      rate2[j*4+0]=rate2[j*4+0]+1.0;
		    }
		  // 2.3.3. sample phi_1
		  U1=runif(0,1);
		  temp=U1*pnorm(phi2[j*I+1]/sdphi21[j], 0, 1, 1, 0)+pnorm(-phi2[j*I+1]/sdphi21[j], 0, 1, 1, 0);
		  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi21[j]+phi2[j*I+1];	
		  
		  Cphi1=pnorm(phi2[j*I+1]/sdphi21[j],0,1,1,1)-pnorm(propphi1/sdphi21[j],0,1,1,1);
		  sumgammaphi1=0;
		  for (s=0;s<S;s++)
		    {
		      if (stateXnrep2[i2*S+s]==1&dataset2[j*S+s]>0)
			{
			  for (j2=0;j2<dataset2[j*S+s];j2++)
			    {
			      sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi2[j*I+1]+j2);
			    }
			}
		    }
		  logAphi1=sumgammaphi1+sumnrep2[i2]*(log(mu2[j*I+1]+phi2[j*I+1])-log(mu2[j*I+1]+propphi1))+nnrep21[i2]*(propphi1*log(propphi1/(mu2[j*I+1]+propphi1))-phi2[j*I+1]*log(phi2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1])))+(Aphi2[j*I+1]-1.0)*log(propphi1/phi2[j*I+1])-Bphi2[j*I+1]*(propphi1-phi2[j*I+1])+Cphi1;
		  U=runif(0,1);
		  if(log(U)<logAphi1)
		    {
		      phi2[j*I+1]=propphi1;
		      rate2[j*4+1]=rate2[j*4+1]+1.0;
		    }
		}
	      // 2.4. sample lambda when method="poisson"(*met=1)
	      if (*met==1)
		{
		  a=Alambda2[j*I+0]+sumZnrep2[i2];
		  b=1/(Blambda2[j*I+0]+nZnrep21[i2]);
		  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda2[0], sumZnrep2[i2], Blambda2[0], nZnrep21[i2]);
		  lambda2[j*I+0]=rgamma(a,b);
		  a=Alambda2[j*I+1]+sumnrep2[i2];
		  b=1/(Blambda2[j*I+1]+nnrep21[i2]);
		  lambda2[j*I+1]=rgamma(a,b);
		  //Rprintf("sample lambda1 %lf %lf %lf %d\n", Alambda2[1], sumnrep2[i2], Blambda2[1], nnrep21[i2]);
		  if(lambda2[j*I+0]>lambda2[j*I+1])
		    {
		      temp=lambda2[j*I+1];
		      lambda2[j*I+1]=lambda2[j*I+0];
		      lambda2[j*I+0]=temp;
		      //Rprintf("lambda= %lf %lf\n", lambda2[0], lambda2[1]);
		    }
		  //Rprintf("lambda= %lf %lf\n", lambda2[0], lambda2[1]);
		}
	      // 2.5. sample mu and phi when method="NB" (*met=2)
	      if (*met==2) // 
		{
		  // 2.5.1.Sample proposal mu from truncated normal distributions (postive tail)
		  U1=runif(0,1);
		  temp=U1*pnorm(mu2[j*I+0]/sdmu20[j], 0, 1, 1, 0)+pnorm(-mu2[j*I+0]/sdmu20[j], 0, 1, 1, 0);
		  propmu0=qnorm(temp, 0, 1, 1, 0)*sdmu20[j]+mu2[j*I+0];	
		  U1=runif(0,1);
		  temp=U1*pnorm(mu2[j*I+1]/sdmu21[j], 0, 1, 1, 0)+pnorm(-mu2[j*I+1]/sdmu21[j], 0, 1, 1, 0);
		  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu21[j]+mu2[j*I+1];	
		  
		  // 2.5.2. Calculate the acceptance ratio for mu
		  Cmu0=pnorm(mu2[j*I+0]/sdmu20[j],0,1,1,1)-pnorm(propmu0/sdmu20[j],0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
		  Cmu1=pnorm(mu2[j*I+1]/sdmu21[j],0,1,1,1)-pnorm(propmu1/sdmu21[j],0,1,1,1);
		  logAmu0=sumZnrep2[i2]*(log(propmu0/(propmu0+phi2[j*I+0]))-log(mu2[j*I+0]/(mu2[j*I+0]+phi2[j*I+0])))+phi2[j*I+0]*nZnrep21[i2]*(log(mu2[j*I+0]+phi2[j*I+0])-log(propmu0+phi2[j*I+0]))+(Amu2[j*I+0]-1.0)*log(propmu0/mu2[j*I+0])-Bmu2[j*I+0]*(propmu0-mu2[j*I+0])+Cmu0;
		  logAmu1=sumnrep2[i2]*(log(propmu1/(propmu1+phi2[j*I+1]))-log(mu2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1])))+phi2[j*I+1]*nnrep21[i2]*(log(mu2[j*I+1]+phi2[j*I+1])-log(propmu1+phi2[j*I+1]))+(Amu2[j*I+1]-1.0)*log(propmu1/mu2[j*I+1])-Bmu2[j*I+1]*(propmu1-mu2[j*I+1])+Cmu1;
		  // 2.5.3. sample mu
		  U=runif(0,1);
		  if(log(U)<logAmu0)
		    {
		      mu2[j*I+0]=propmu0;
		      rate2[j*4+2]=rate2[j*4+2]+1.0;
		    }
		  U=runif(0,1);
		  if(log(U)<logAmu1)
		    {
		      mu2[j*I+1]=propmu1;
		      rate2[j*4+0]=rate2[j*4+0]+1.0;
		    }
		  // 2.5.4. Sample proposal phi from truncated normal distribution (postive tail)
		  U1=runif(0,1);
		  temp=U1*pnorm(phi2[j*I+0]/sdphi20[j], 0, 1, 1, 0)+pnorm(-phi2[j*I+0]/sdphi20[j], 0, 1, 1, 0);
		  propphi0=qnorm(temp, 0, 1, 1, 0)*sdphi20[j]+phi2[j*I+0];	
		  U1=runif(0,1);
		  temp=U1*pnorm(phi2[j*I+1]/sdphi21[j], 0, 1, 1, 0)+pnorm(-phi2[j*I+1]/sdphi21[j], 0, 1, 1, 0);
		  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi21[j]+phi2[j*I+1];	
		  
		  // 2.5.5. Calculate the acceptance ratio for phi
		  Cphi0=pnorm(phi2[j*I+0]/sdphi20[j],0,1,1,1)-pnorm(propphi0/sdphi20[j],0,1,1,1);// log(PHI(phi^t/sd_phi))-log(PHI(phi'/sd_phi))##using R function pnorm(x, mu, phi, iflowtail=1, iflog=1)
		  Cphi1=pnorm(phi2[j*I+1]/sdphi21[j],0,1,1,1)-pnorm(propphi1/sdphi21[j],0,1,1,1);
		  sumgammaphi0=0;
		  sumgammaphi1=0;
		  for (s=0;s<S;s++)
		    {
		      if (stateXnrep2[i2*S+s]==0&stateZ2[j*S+s]==1&dataset2[j*S+s]>0)
			{
			  for (j2=0;j2<dataset2[j*S+s];j2++)
			    {
			      sumgammaphi0=sumgammaphi0+log(propphi0+j2)-log(phi2[j*I+0]+j2);
			    }
			}
		      if (stateXnrep2[i2*S+s]==1&dataset2[j*S+s]>0)
			{
			  for (j2=0;j2<dataset2[j*S+s];j2++)
			    {
			      sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi2[j*I+1]+j2);
			    }
			}	
		    }
		  logAphi0=sumgammaphi0+sumZnrep2[i2]*(log(mu2[j*I+0]+phi2[j*I+0])-log(mu2[j*I+0]+propphi0))+nZnrep21[i2]*(propphi0*log(propphi0/(mu2[j*I+0]+propphi0))-phi2[j*I+0]*log(phi2[j*I+0]/(mu2[j*I+0]+phi2[j*I+0])))+(Aphi2[j*I+0]-1.0)*log(propphi0/phi2[j*I+0])-Bphi2[j*I+0]*(propphi0-phi2[j*I+0])+Cphi0;
		  logAphi1=sumgammaphi1+sumnrep2[i2]*(log(mu2[j*I+1]+phi2[j*I+1])-log(mu2[j*I+1]+propphi1))+nnrep21[i2]*(propphi1*log(propphi1/(mu2[j*I+1]+propphi1))-phi2[j*I+1]*log(phi2[j*I+1]/(mu2[j*I+1]+phi2[j*I+1])))+(Aphi2[j*I+1]-1.0)*log(propphi1/phi2[j*I+1])-Bphi2[j*I+1]*(propphi1-phi2[j*I+1])+Cphi1;
		  // 2.5.6. sample phi
		  U=runif(0,1);
		  if(log(U)<logAphi0)
		    {
		      phi2[j*I+0]=propphi0;
		      rate2[j*4+3]=rate2[j*4+3]+1.0;
		    }
		  U=runif(0,1);
		  if(log(U)<logAphi1)
		    {
		      phi2[j*I+1]=propphi1;
		      rate2[j*4+1]=rate2[j*4+1]+1.0;
		    }
		  //Rprintf("mu,phi= %lf %lf %lf %lf\n",mu2[j*I+1], phi2[j*I+1], mu2[j*I+0],  phi2[j*I+0]);
		  if(mu2[j*I+0]>mu2[j*I+1])
		    {
		      //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu2[0], mu2[1], phi2[0], phi2[1]);
		      temp=phi2[j*I+1];
		      phi2[j*I+1]=phi2[j*I+0];
		      phi2[j*I+0]=temp;
		      temp=mu2[j*I+1];
		      mu2[j*I+1]=mu2[j*I+0];
		      mu2[j*I+0]=temp;
		    }
		  //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu2[0], mu2[1], phi2[0], phi2[1]);
		}
	    }
	}
	
      // 0.1.0.2. Calculate the acceptance ratio for c
      Cc=pnorm(c/sdc,0,1,1,1)-pnorm(propc/sdc,0,1,1,1);//formula is log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu)), where PHI is standard norm cdf##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
      logAc=logAc1-logAc2+Cc;
      // 0.1.0.3. sample c
      U=runif(0,1);
      if(log(U)<logAc)
	{
	  c=propc;
	  rate[nrate-1]=rate[nrate-1]+1.0;
	}
      

      // iteration step 3. calculate likelihood function 
      logl=0.0;  
      if(nsr1>0)
	{
	  pX1[0]=1/(1+exp(log(1-qrep1[1])-log(qrep1[0])));//stationary initial distribution
	  pX0[0]=1-pX1[0];
	  for (s=0; s<S; s++)
	    {
	      for (j=0;j<nsr1;j++)
		{
		  logl=logl+log(((1-pi1[j])*indexY1[j*S+s]+pi1[j]*pY1[0][j*S+s]*indexY1[j*S+s]+pY1[0][j*S+s]*(1-indexY1[j*S+s]))*pX1[0]+pY1[1][j*S+s]*pX1[1]);
		}
	    }
	}
      if (nsr2>0)
	{
	  pX1[nsp1+indexnsr1+0]=1/(1+exp(log(1-qrep2[1])-log(qrep2[0])));//stationary initial distribution
	  pX0[nsp1+indexnsr1+0]=1-pX1[nsp1+indexnsr1+0];
	  for (s=0; s<S; s++)
	    {
	      for (j=0;j<nsr2;j++)
		{
		  logl=logl+log(((1-pi2[j])*indexY2[j*S+s]+pi2[j]*pY2[0][j*S+s]*indexY2[j*S+s]+pY2[0][j*S+s]*(1-indexY2[j*S+s]))*pX1[nsp1+indexnsr1+0]+pY2[1][j*S+s]*pX1[nsp1+indexnsr1+1]);
		}
	    }
	}
      if (nsp1>0)
	{
	  for (j=nsr1;j<nr1;j++)
	    {     
	      i1=j-nsr1;
	      pX1[i1+indexnsr1]=1/(1+exp(log(1-qnrep1[i1*I+1])-log(qnrep1[i1*I+0])));//stationary initial distribution
	      pX0[i1+indexnsr1]=1-pX1[i1+indexnsr1];
	    }
	  for (s=0; s<S; s++)
	    {
	      for (j=nsr1;j<nr1;j++)
		{ 
		  i1=j-nsr1;
		  logl=logl+log(((1-pi1[j])*indexY1[j*S+s]+pi1[j]*pY1[0][j*S+s]*indexY1[j*S+s]+pY1[0][j*S+s]*(1-indexY1[j*S+s]))*pX0[i1+indexnsr1]+pY1[1][j*S+s]*pX1[i1+indexnsr1]);
		}
	    }
	}
      if (nsp2>0)
	{
	  for (j=nsr2;j<nr2;j++)
	    {     
	      i2=j-nsr2;
	      pX1[nsp1+indexnsr1+indexnsr2+i2]=1/(1+exp(log(1-qnrep2[i2*I+1])-log(qnrep2[i2*I+0])));//stationary initial distribution
	      pX0[nsp1+indexnsr1+indexnsr2+i2]=1-pX1[nsp1+indexnsr1+indexnsr2+i2];
	    }
	  for (s=0; s<S; s++)
	    {
	      for (j=nsr2;j<nr2;j++)
		{ 
		  i2=j-nsr2;
		  logl=logl+log(((1-pi2[j])*indexY2[j*S+s]+pi2[j]*pY2[0][j*S+s]*indexY2[j*S+s]+pY2[0][j*S+s]*(1-indexY2[j*S+s]))*pX0[nsp1+indexnsr1+indexnsr2+i2]+pY2[1][j*S+s]*pX1[nsp1+indexnsr1+indexnsr2+i2]);
		}
	    }
	}
      //Rprintf("loglikelihood= %lf\n", logl);

      // iteration step 4. pass the results back to R
      if(nsr1>0)
	{
	  for (j=0;j<nsr1;j++)
	    {
	      q1[j*I+0]=qrep1[0];
	      q1[j*I+1]=qrep1[1];
	    }
	}
      if (nsp1>0)
	{
	  for (j=0;j<nsp1;j++)
	    {
	      q1[(j+nsr1)*I+0]=qnrep1[j*I+0];
	      q1[(j+nsr1)*I+1]=qnrep1[j*I+1];
	    }
	}
      if (nsr2>0)
	{
	  for (j=0;j<nsr2;j++)
	    {
	      q2[j*I+0]=qrep2[0];
	      q2[j*I+1]=qrep2[1];
	    }
	}
      if (nsp2>0)
	{
	  for (j=0;j<nsp2;j++)
	    {
	      q2[(j+nsr2)*I+0]=qnrep2[j*I+0];
	      q2[(j+nsr2)*I+1]=qnrep2[j*I+1];
	    }
	}
      
      if (z>=bN&z==j1*10+bN)
	{
	  for (j=0;j<nr1;j++)
	    {
	      *(es_q11+j1*nr1+j)=q1[j*I+1];
	      *(es_q10+j1*nr1+j)=q1[j*I+0];
	      *(es_pi1+j1*nr1+j)=pi1[j];
	      *(es_lambda11+j1*nr1+j)=lambda1[j*I+1];
	      *(es_lambda10+j1*nr1+j)=lambda1[j*I+0];
	      *(es_mu11+j1*nr1+j)=mu1[j*I+1];
	      *(es_phi11+j1*nr1+j)=phi1[j*I+1];
	      *(es_mu10+j1*nr1+j)=mu1[j*I+0];
	      *(es_phi10+j1*nr1+j)=phi1[j*I+0];
	    }
	  for (j=0;j<nr2;j++)
	    {
	      *(es_q21+j1*nr2+j)=q2[j*I+1];
	      *(es_q20+j1*nr2+j)=q2[j*I+0];
	      *(es_pi2+j1*nr2+j)=pi2[j];
	      *(es_lambda21+j1*nr2+j)=lambda2[j*I+1];
	      *(es_lambda20+j1*nr2+j)=lambda2[j*I+0];
	      *(es_mu21+j1*nr2+j)=mu2[j*I+1];
	      *(es_phi21+j1*nr2+j)=phi2[j*I+1];
	      *(es_mu20+j1*nr2+j)=mu2[j*I+0];
	      *(es_phi20+j1*nr2+j)=phi2[j*I+0];
	    }
	  *(loglikeli+j1)=logl;  
	  j1++;
	}
      
      if (z==z1*jN)
	{
	  Rprintf("%d % ",z1*10);
	  z1++;	    
	}
      //if (z==j1*1000)
      //	{
      //j1=j1++;
      /*
          Rprintf("MCMC step= %d \n", z+1);
	  for (j=0;j<nr1;j++)
	    {
	      Rprintf("parameters of %d", j+1);
	      Rprintf("th experiment for first condition are\n");
	      Rprintf("%lf %lf %lf %lf %lf \n", mu1[j*I+1], phi1[j*I+1], pi1[j], mu1[j*I+0], phi1[j*I+0]);
	      Rprintf("and the transition probability and ratio= %lf %lf %lf %lf\n", q1[j*I+1], q1[j*I+0],c, q1[j*I+0]/(q1[j*I+0]+1-q1[j*I+1]));
	    }
	  for (j=0;j<nr2;j++)
	    {
	      Rprintf("parameters of %d", j+1);
	      Rprintf("th experiment for second condition are \n");
	      Rprintf(" %lf %lf %lf %lf %lf\n", mu2[j*I+1], phi2[j*I+1], pi2[j], mu2[j*I+0], phi2[j*I+0]);
	      Rprintf("and the transition probability and ratio= %lf %lf %lf %lf\n", q2[j*I+1], q2[j*I+0],c, q2[j*I+0]/(q2[j*I+0]+1-q2[j*I+1]));
	    }
	    //	}*/
    }
  Rprintf("\n");
  PutRNGstate();
  double iter;
  iter=*N;  
  for (i=0;i<nr1; i++)
    {
      *(acrate1+i*4+0)=rate1[i*4+0]/iter;
      *(acrate1+i*4+1)=rate1[i*4+1]/iter;
      *(acrate1+i*4+2)=rate1[i*4+2]/iter;
      *(acrate1+i*4+3)=rate1[i*4+3]/iter;
      //Rprintf("The acceptance rate for %d", i+1);
      //Rprintf("th experiment for first condition are:\n");
      //Rprintf("%lf %lf %lf %lf \n",rate1[i*4+0]/iter,rate1[i*4+1]/iter, rate1[i*4+2]/iter, rate1[i*4+3]/iter);// acceptance rate
    }
  for (i=0;i<nr2; i++)
    {
      *(acrate2+i*4+0)=rate2[i*4+0]/iter;
      *(acrate2+i*4+1)=rate2[i*4+1]/iter;
      *(acrate2+i*4+2)=rate2[i*4+2]/iter;
      *(acrate2+i*4+3)=rate2[i*4+3]/iter;
      //Rprintf("The acceptance rate for %d", i+1);
      //Rprintf("th experiment for second condition are:\n");
      //Rprintf("%lf %lf %lf %lf \n",rate2[i*4+0]/iter,rate2[i*4+1]/iter, rate2[i*4+2]/iter, rate2[i*4+3]/iter);// acceptance rate
    }
  
  //Rprintf("The acceptance rate for transition probabilities are:\n");
  for (i=0;i<nrate;i++)
    { 
      *(acrate+i)=rate[i]/iter;
      //Rprintf("%lf ",rate[i]/iter);// acceptance rate
    }
  //Rprintf("\n");
 
  if (nsr1>0)
    {
      for (i=0;i<nsr1;i++)
	{
	  for (s=0; s<S; s++)
	    {
	      *(PP1+i*S+s)=sumprobXrep1[s]/iter;
	    }
	}
    }
  if (nsp1>0)
    {
      for (i=nsr1;i<nr1;i++)
	{
	  i1=i-nsr1;
	  for (s=0;s<S;s++)
	    {
	      *(PP1+i*S+s)=sumprobXnrep1[i1*S+s]/iter;
	    }
	}
    }
  if (nsr2>0)
    {
      for (i=0;i<nsr2;i++)
	{
	  for (s=0; s<S; s++)
	    {
	      *(PP2+i*S+s)=sumprobXrep2[s]/iter;
	    }
	}
    }
  if (nsp2>0)
    {
      for (i=nsr2;i<nr2;i++)
	{
	  i2=i-nsr2;
	  for (s=0;s<S;s++)
	    {
	      *(PP2+i*S+s)=sumprobXnrep2[i2*S+s]/iter;
	    }
	}
    }
   for (i=0; i<I; i++)
    {
      free(pY1[i]);
      free(logpY1[i]);
      free(pY2[i]);
      free(logpY2[i]);
      if (nsr1>0)
	{
	  free(pXcYrep1[i]);
	}
      if (nsr2>0)
	{
	  free(pXcYrep2[i]);
	}
      if (nsp1>0)
	{
	  free(pXcYnrep1[i]);
	}
      if (nsp2>0)
	{
	  free(pXcYnrep2[i]);
	}
    }
  free(pY1);
  free(logpY1);
  free(pY2);
  free(logpY2);
  free(stateZ1);
  free(stateZ2);
  if (nsr1>0)
    {
      free(pXcYrep1);
      free(probXrep1);
      free(sumprobXrep1);
      free(stateXrep1);
    }
  if (nsr2>0)
    {
      free(pXcYrep2);
      free(probXrep2);
      free(sumprobXrep2);
      free(stateXrep2);
    }
  if (nsp1>0)
    {
      free(pXcYnrep1);
      free(probXnrep1);
      free(sumprobXnrep1);
      free(stateXnrep1);
    }
  if (nsp2>0)
    {
      free(pXcYnrep2);
      free(probXnrep2);
      free(sumprobXnrep2);
      free(stateXnrep2);
    }
}
