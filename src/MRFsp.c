//***********************************************************//
//***********************************************************//
//***    The program is used to do ChIP-Seq data analysis ***//
//***    for histone modification using Markov Random     ***//
//***    Field (MRF model), which will form the core part ***//
//***    of package enRich. The program is written by      ***//
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
void MRFsp(int *data1, int *nsp, int *size, int *met, double *qprior, double *piprior, double *poisprior, double *NBprior, int *N, int *Nb, int *jumpN, double *var, double *var1, double *PP, double *es_pi, double *es_q1, double *es_q0, double *es_lambda1, double *es_mu1, double *es_phi1, double *es_lambda0, double *es_mu0, double *es_phi0, double *loglikeli, double *acrate, double *acrate1)
{
  int S=*size;
//  int N1=*N;
  int bN=*Nb;
  int jN=*jumpN;
  int nr=*nsp;
  int nrI=nr*I;
  int nrS=nr*S;
  //Rprintf("S, N1, nr, nrI and nrS= %d %d %d %d %d\n", S, N1, nr, nrI, nrS);
  double *AQ, *BQ;  // priors of transition probabilities q0.
  AQ=(double*)malloc(sizeof(double*)*nr);
  BQ=(double*)malloc(sizeof(double*)*nr);
  double *Api, *Bpi;// priors of mixture portion of ZIP model.
  Api=(double*)malloc(sizeof(double*)*nr);
  Bpi=(double*)malloc(sizeof(double*)*nr);
  double *Alambda, *Blambda; // priors of lambda for poisson fitting.
  Alambda=(double*)malloc(sizeof(double*)*nrI);
  Blambda=(double*)malloc(sizeof(double*)*nrI);
  double *Amu, *Bmu;//priors of mu_i, which is the mean of NB for ith component;
  Amu=(double*)malloc(sizeof(double*)*nrI);
  Bmu=(double*)malloc(sizeof(double*)*nrI);
  double *Aphi, *Bphi;//priors of phi_i, which is the overdispersion of NB for ithe component;
  Aphi=(double*)malloc(sizeof(double*)*nrI);
  Bphi=(double*)malloc(sizeof(double*)*nrI);
  int *data;// bin counts--observations.
  data=(int*)malloc(sizeof(int*)*nrS);
  int *stateX; // bin states, latent variables
  stateX=(int*)malloc(sizeof(int*)*nrS);
  int *stateZ; // inner latent variables
  stateZ=(int*)malloc(sizeof(int*)*nrS);
  double *lambda, *mu, *phi, *q, *pi; // parameters of Poisson-Gamma mixing distribution and transition probabilities
  lambda=(double*)malloc(sizeof(double*)*nrI);
  mu=(double*)malloc(sizeof(double*)*nrI);
  phi=(double*)malloc(sizeof(double*)*nrI);
  q=(double*)malloc(sizeof(double*)*nrI);
  pi=(double*)malloc(sizeof(double*)*nr);
  int i, j, j1, j2, s, z; 
 
  double *probX, *sumprobX, *probZ, *pX1, *pX0;
  int *indexY;
  probX=(double*)malloc(sizeof(double*)*nrS);
  sumprobX=(double*)malloc(sizeof(double*)*nrS);
  probZ=(double*)malloc(sizeof(double*)*nr); 
  pX1=(double*)malloc(sizeof(double*)*nr); 
  pX0=(double*)malloc(sizeof(double*)*nr);
  indexY=(int*)malloc(sizeof(int*)*nrS);

  int n[I][I];//the counts of pair
  double Q[I][I],**logpY,**pY,**pXcY; 

  logpY=(double**)malloc(sizeof(double*)*I);
  pY=(double**)malloc(sizeof(double*)*I);
  pXcY=(double**)malloc(sizeof(double*)*I);
 
  for(i=0;i<I;i++)
    {
      logpY[i]=(double*)malloc(sizeof(double*)*nrS);
      pY[i]=(double*)malloc(sizeof(double*)*nrS);
      pXcY[i]=(double*)malloc(sizeof(double*)*nrS);
    }
  
  double *sum2X;//using for counting pairs of X
  sum2X=(double*)malloc(sizeof(double*)*(S-1));
  double U,U1,temp,a,b,logl;
  
  ///////// R data1 pass to data
  for(s=0;s<nrS;s++)
    {
      data[s]=*(data1+s);
      sumprobX[s]=0.0;
    }
  //  for(s=47;s<55;s++)//48644;s<48650;s++)
  //{
  //Rprintf("data[48644:48650,]= %d %d\n", data[s], data[S+s]);
  //}

  ///////// give priors for parameters; 
 
  for (i=0;i<nr;i++)
    {
      AQ[i]=*(qprior+i*I+0);
      BQ[i]=*(qprior+i*I+1);
    
      Api[i]=*(piprior+i*I+0);
      Bpi[i]=*(piprior+i*I+1);
  
      Alambda[i*I+1]=*(poisprior+i*4+0);
      Blambda[i*I+1]=*(poisprior+i*4+1);
      Alambda[i*I+0]=*(poisprior+i*4+2);
      Blambda[i*I+0]=*(poisprior+i*4+3);

      Amu[i*I+1]=*(NBprior+i*8+0);
      Bmu[i*I+1]=*(NBprior+i*8+1);
      Aphi[i*I+1]=*(NBprior+i*8+2);
      Bphi[i*I+1]=*(NBprior+i*8+3);
      Amu[i*I+0]=*(NBprior+i*8+4);
      Bmu[i*I+0]=*(NBprior+i*8+5);
      Aphi[i*I+0]=*(NBprior+i*8+6);
      Bphi[i*I+0]=*(NBprior+i*8+7);

      //Rprintf("prior of q= %f %f %f %f \n",AQ[i*I+0], AQ[i*I+1], BQ[i*I+0], BQ[i*I+1]);
      //Rprintf("Prior of %d",i);
      //Rprintf("th experiment, prior of pi= %f %f \n", Api[i], Bpi[i]);
      //Rprintf("prior of lambda= %f %f %f %f \n", Alambda[i*I+1], Blambda[i*I+1], Alambda[i*I+0], Blambda[i*I+0]);
      //Rprintf("prior of mu and phi= %f %f %f %f %f %f %f %f \n", Amu[i*I+1], Bmu[i*I+1], Aphi[i*I+1], Bphi[i*I+1], Amu[i*I+0], Bmu[i*I+0], Aphi[i*I+0], Bphi[i*I+0]); 
    }
  
  //////// give initial values, set X=1 if Y>cp, where cp is initial cutting points
  double *sumZ1, *sum1;
  sumZ1=(double*)malloc(sizeof(double*)*nr);
  sum1=(double*)malloc(sizeof(double*)*nr);
  int *n1,*n0, *nZ1,*nZ0;
  n1=(int*)malloc(sizeof(int*)*nr);
  n0=(int*)malloc(sizeof(int*)*nr);
  nZ1=(int*)malloc(sizeof(int*)*nr);
  nZ0=(int*)malloc(sizeof(int*)*nr);

  double c=1; //the ratio of transtion probability
  for (i=0;i<nr;i++)
    {   
      sumZ1[i]=0.0;
      sum1[i]=0.0;
      n1[i]=0;
      n0[i]=0;
      nZ1[i]=0;
      nZ0[i]=0;
      for (s=0;s<S;s++)
	{ 
	   if(data[i*S+s]>cp)
	    { 
	      stateX[i*S+s]=1;
	      n1[i]++;
	      sum1[i]=sum1[i]+data[i*S+s];
	      stateZ[i*S+s]=2;
	      indexY[i*S+s]=0;
	    }
	  else
	    {
	      stateX[i*S+s]=0;
	      n0[i]++;
	      if (data[i*S+s]==0)
		{
		  indexY[i*S+s]=1;
		  stateZ[i*S+s]=0;
		  nZ0[i]++;
		}
	      else
		{
		  indexY[i*S+s]=0;
		  stateZ[i*S+s]=1;
		  sumZ1[i]=sumZ1[i]+data[i*S+s];
		  nZ1[i]++;
		}
	    }
	}
      if ((sum1[i]==0)|(n1[i]==0))
	{
	  sum1[i]=5;
	  n1[i]=1;
	}
      //printf("initial sumZ1, sum1, sum,nZ0, nZ1,nZ0+nZ1, n1, n0, n1+n0=% lf %lf %d %d %d %d %d %d \n",sumZ1[i], sum1[i], nZ0[i], nZ1[i], nZ0[i]+nZ1[i], n0[i], n1[i], n0[i]+n1[i]);
      if (*met==0)
	{
	  lambda[i*I+0]=sumZ1[i]/nZ1[i];
	  lambda[i*I+1]=0;
	  mu[i*I+1]=sum1[i]/n1[i];
	  mu[i*I+0]=0;
	  phi[i*I+1]=1;
	  phi[i*I+0]=0;
	}
      if (*met==1)
	{
	  lambda[i*I+0]=sumZ1[i]/nZ1[i];//sum0/n0;
	  lambda[i*I+1]=sum1[i]/n1[i]; 
	  mu[i*I+0]=0;
	  mu[i*I+1]=0;
	  phi[i*I+0]=0;
	  phi[i*I+1]=0;
	  //printf("initial lambda %lf %lf\n", lambda[0],lambda[1]); 
	}
      if (*met==2)
	{
	  lambda[i*I+0]=0;//sum0/n0;
	  lambda[i*I+1]=0; 
	  mu[i*I+0]=sumZ1[i]/nZ1[i];
	  mu[i*I+1]=sum1[i]/n1[i];
	  phi[i*I+0]=1;
	  phi[i*I+1]=1;
	  //printf("initial mu and phi %lf %lf %lf %lf\n", mu[i*I+1], phi[i*I+1],mu[i*I+0], phi[i*I+0]);
	}
      pi[i]=(nZ1[i]+0.0)/(nZ0[i]+0.0);
      //printf("initial nZ and pi %d %d %lf\n", nZ1[i], nZ0[i], pi[i]);
      q[i*I+0]=0.1; //prob(X_s=1|X_{s-1}=0)
      q[i*I+1]=0.9; //prob(X_s=1|X_{s-1}=1)
    }
  
  //for(s=48605;s<48610;s++)
  //{
  //  Rprintf("indexY[48606:48610,]= %d %d\n", indexY[s], indexY[S+s]);
  //}
  
  //Rprintf("initial values\n");
  //for (j=0;j<nr;j++)
  //{
  //  Rprintf("parameters of %d", j+1);
  //  Rprintf("th experiment are %lf %lf %lf %lf %lf %lf %lf\n", q[j*I+0], q[j*I+1],  mu[j*I+1], phi[j*I+1],pi[j], mu[j*I+0], phi[j*I+0]);
  //}  

  double *sdmu0, *sdmu1, *sdphi0, *sdphi1, *sdq0, sdc, *rate, *rate1, sumgammaphi0, sumgammaphi1;
  sdmu0=(double*)malloc(sizeof(double*)*nr);
  sdmu1=(double*)malloc(sizeof(double*)*nr);
  sdphi0=(double*)malloc(sizeof(double*)*nr); 
  sdphi1=(double*)malloc(sizeof(double*)*nr); 
  sdq0=(double*)malloc(sizeof(double*)*nr);
  int nr4=4*nr;
  rate=(double*)malloc(sizeof(double*)*nr4);
  int nrplus1=nr+1;
  rate1=(double*)malloc(sizeof(double*)*nrplus1);
  for (i=0;i<nr;i++)
    {
      sdmu1[i]=*(var+i*4+0);
      sdphi1[i]=*(var+i*4+1);
      sdmu0[i]=*(var+i*4+2);
      sdphi0[i]=*(var+i*4+3);
      sdq0[i]=*(var1+i);
      rate[i*4+0]=0.0;
      rate[i*4+1]=0.0;
      rate[i*4+2]=0.0;
      rate[i*4+3]=0.0;
      rate1[i]=0.0;
      //Rprintf("var= %lf %lf %lf %lf \n", sdmu1[i], sdphi1[i], sdmu0[i], sdphi0[i]);
    }
  rate1[nr]=0.0;
  sdc=*(var1+nr);
  double propmu0, propmu1, propphi0, propphi1, Cmu0, Cmu1, logAmu0, logAmu1, Cphi0, Cphi1, logAphi0, logAphi1;
  double propc, Cc, logAc, logAc1, logAc2, propq0, Cq0, logAq0;

  j1=0;
  int z1=0;
  GetRNGstate();
  for(z=0;z<*N;z++)//show results for every 10th iterations
    {
      // 2.1.0 sample the common protion c for transition probabilities. We sample c by Metropolis-Hasting method
      // 2.1.0.1 sample proposal value for c
      U1=runif(0,1);
      temp=U1*pnorm(c/sdc, 0, 1, 1, 0)+pnorm(-c/sdc, 0, 1, 1, 0);
      propc=qnorm(temp, 0, 1, 1, 0)*sdc+c;
      logAc1=0;
      logAc2=0;
      
      for (j=0;j<nr;j++)
	{	    
	  n[0][0]=0;
	  n[0][1]=0;
	  n[1][0]=0;
	  n[1][1]=0;
	  Q[0][0]=1-q[j*I+0];
	  Q[0][1]=q[j*I+0];
	  Q[1][0]=1-q[j*I+1];
	  Q[1][1]=q[j*I+1];	 
	  n0[j]=0;
	  n1[j]=0;
	  sumZ1[j]=0.0;
	  sum1[j]=0.0;
	  nZ0[j]=0;
	  nZ1[j]=0;
	  // iteration step 1: Given parameters renew the states X of bins by using direct gibbs method. We also sample the inner states Z given X==0. 
	  // 1.1 sample stateX
	  for (s=0; s<S; s++)
	    {
	      if (*met==0)// poissonNB
		{
		  logpY[0][j*S+s]=dpois(data[j*S+s], lambda[j*I+0], 1);//use R function dpois(y, lambda, logistrue=1);
		  pY[0][j*S+s]=dpois(data[j*S+s], lambda[j*I+0], 0);//use R function dpois(y, lambda, logisfalse=0);
		  logpY[1][j*S+s]=dnbinom(data[j*S+s], phi[j*I+1], phi[j*I+1]/(mu[j*I+1]+phi[j*I+1]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
		  pY[1][j*S+s]=dnbinom(data[j*S+s], phi[j*I+1], phi[j*I+1]/(mu[j*I+1]+phi[j*I+1]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
		}
	      if (*met==1)// poisson
		{
		  for (i=0;i<I;i++)
		    {
		      logpY[i][j*S+s]=dpois(data[j*S+s], lambda[j*I+i], 1);//use R function dpois(y, lambda, logistrue=1);
		      pY[i][j*S+s]=dpois(data[j*S+s], lambda[j*I+i], 0);//use R function dpois(y, lambda, logisfalse=0);
		    }
		}
	      if (*met==2)// NB
		{
		  for (i=0;i<I;i++)
		    {
		      logpY[i][j*S+s]=dnbinom(data[j*S+s], phi[j*I+i], phi[j*I+i]/(mu[j*I+i]+phi[j*I+i]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
		      pY[i][j*S+s]=dnbinom(data[j*S+s], phi[j*I+i], phi[j*I+i]/(mu[j*I+i]+phi[j*I+i]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
		    }
		}
	    }
	
	  //Rprintf("mu and phi are: %lf %lf %lf %lf %lf %lf %lf %lf\n", mu[1], phi[1], mu[0], phi[0], mu[3], phi[3], mu[2], mu[2]);
	  //for(s=47;s<51;s++)
	  //{
	  //Rprintf("pY[47:51,]= %lf %lf %lf %lf \n", pY[1][s], pY[0][s], pY[1][S+s], pY[0][S+s]);
	  //Rprintf("Y and sumprobX[47:55,]= %d %d %lf %lf \n", data[s], data[S+s], sumprobX[s], sumprobX[S+s]);
	  //}
	  pXcY[0][j*S+0]=((1-pi[j])*indexY[j*S+0]+pi[j]*indexY[j*S+0]*pY[0][j*S+0]+(1-indexY[j*S+0])*pY[0][j*S+0])*Q[0][stateX[j*S+1]];
	  pXcY[1][j*S+0]=pY[1][j*S+0]*Q[1][stateX[j*S+1]];
	  if ((pXcY[0][j*S+0]<=0)&&(pXcY[1][j*S+0]<=0)&&(data[j*S+0]>cp))
	{ 
		probX[j*S+0]=1;
	}
	else if (pXcY[1][j*S+0]<=0)
	{
		probX[j*S+0]=0;
	}
	else if (pXcY[0][j*S+0]<=0)
	{
		probX[j*S+0]=1;
	}
	else
	{
      	probX[j*S+0]=1/(1+exp(log(pXcY[0][j*S+0])-log(pXcY[1][j*S+0])));//posterior distribution of P(X_1=1)
	}
	  sumprobX[j*S+0]=sumprobX[j*S+0]+probX[j*S+0];
	  U=runif(0,1);      
	  if(U<probX[j*S+0])
	    {
	      stateX[j*S+0]=1;
	      n1[j]++;//To get sum_s I(X_s=1)
	    }
	  else
	    {
	      stateX[j*S+0]=0;
	      n0[j]++;
	    }
	  //sumprobX[j*S+0]=sumprobX[j*S+0]+stateX[j*S+0];
	  for (s=1; s<S-1; s++)
	    {
	      pXcY[0][j*S+s]=((1-pi[j])*indexY[j*S+s]+pi[j]*indexY[j*S+s]*pY[0][j*S+s]+(1-indexY[j*S+s])*pY[0][j*S+s])*Q[stateX[j*S+s-1]][0]*Q[0][stateX[j*S+s+1]];
	      pXcY[1][j*S+s]=pY[1][j*S+s]*Q[stateX[j*S+s-1]][1]*Q[1][stateX[j*S+s+1]];
		if ((pXcY[0][j*S+s]<=0)&&(pXcY[1][j*S+s]<=0)&&(data[j*S+s]>cp))
	  { 
		probX[j*S+s]=1;
	  }
	  else if (pXcY[1][j*S+s]<=0)
	  {
		probX[j*S+s]=0;
	  }
	  else if (pXcY[0][j*S+s]<=0)
	  {
		probX[j*S+s]=1;
	  }
	  else
	  {
	      probX[j*S+s]=1/(1+exp(log(pXcY[0][j*S+s])-log(pXcY[1][j*S+s])));//posterior distribution of P(X_1=1)
	  }
	      sumprobX[j*S+s]=sumprobX[j*S+s]+probX[j*S+s];
	      U=runif(0,1);
	      if(U<probX[j*S+s])
		{
		  stateX[j*S+s]=1;
		  n1[j]++;//To get sum_s I(X_s=1)
		}     
	      else
		{
		  stateX[j*S+s]=0; 
		  n0[j]++;
		}
	      //sumprobX[j*S+s]=sumprobX[j*S+s]+stateX[j*S+s];
	    }
	  pXcY[0][j*S+S-1]=((1-pi[j])*indexY[j*S+S-1]+pi[j]*indexY[j*S+S-1]*pY[0][j*S+S-1]+(1-indexY[j*S+S-1])*pY[0][j*S+S-1])*Q[stateX[j*S+S-2]][0];
	  pXcY[1][j*S+S-1]=pY[1][j*S+S-1]*Q[stateX[j*S+S-2]][1];
	  if ((pXcY[0][j*S+S-1]<=0)&&(pXcY[1][j*S+S-1]<=0)&&(data[j*S+S-1]>cp))
	{ 
		probX[j*S+S-1]=1;
	}
	else if (pXcY[1][j*S+S-1]<=0)
	{
		probX[j*S+S-1]=0;
	}
	else if (pXcY[0][j*S+S-1]<=0)
	{
		probX[j*S+S-1]=1;
	}
	else
	{
	  	probX[j*S+S-1]=1/(1+exp(log(pXcY[0][j*S+S-1])-log(pXcY[1][j*S+S-1])));//posterior distribution of P(X_1=1)
	}
	  sumprobX[j*S+S-1]=sumprobX[j*S+S-1]+probX[j*S+S-1];
	  U=runif(0,1);
	  if(U<probX[j*S+S-1])
	    {
	      stateX[j*S+S-1]=1;
	      n1[j]++;//To get sum_s I(X_s=1)
	    }     
	  else
	    {	  
	      stateX[j*S+S-1]=0;
	      n0[j]++;
	    }
      	  //sumprobX[j*S+S-1]=sumprobX[j*S+S-1]+stateX[j*S+S-1];
	  // 1.2 given stateX we sample stateZ
    	  if (*met<2)
	    {
	      probZ[j]=1/(1+exp(log(1-pi[j])-log(pi[j])+lambda[j*I+0]));
	    }
	  if (*met==2)
	    {
	      probZ[j]=1/(1+exp(log(1-pi[j])-log(pi[j])-phi[j*I+0]*log(phi[j*I+0]/(phi[j*I+0]+mu[j*I+0]))));
	    }
	  //Rprintf("probZ= %lf\n", probZ[j]);
	  for (s=0;s<S;s++)
	    {
	      if (stateX[j*S+s]==0)
		{
		  if (data[j*S+s]>0)
		    {
		      stateZ[j*S+s]=1;
		      sumZ1[j]=sumZ1[j]+data[j*S+s];
		      nZ1[j]++;
		    }
		  else
		    {
		      U=runif(0,1);
		      if (U<probZ[j])
			{
			  stateZ[j*S+s]=1;
			  sumZ1[j]=sumZ1[j]+data[j*S+s]; //To get sum_s Y_sI(X_s=0, Z_s=1)
			  nZ1[j]++;//To get sum_s I(X_s=0, Z_s=1)
			}
		      else
			{
			  stateZ[j*S+s]=0;
			  nZ0[j]++;//To get sum_s I(X_s=0, Z_s=0)
			}
		    }
		}
	      else
		{		  
		  sum1[j]=sum1[j]+data[j*S+s]; //To get sum_s Y_sI(X_s=1)
		  stateZ[j*S+s]=2;
		}
	    }
	
          // iteration step 2: sample the parameters from posterior distribution.
	  // 2.1 sample transition probabilities from posterior beta distribution (this will be done after all experiments)
	  // 2.1.1 count the pair
	  for (s=0;s<S-1;s++)
	    {
	      sum2X[s]=stateX[j*S+s]+stateX[j*S+s+1]*2;
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
	  
	  // 2.1.2 sample transition probability q_1, and later sample the common protion c. We sample q_1 by Metropolis-Hasting method
	  U1=runif(0,1);
	  temp=U1*(pnorm((1-q[j*I+0])/sdq0[j], 0, 1, 1, 0)-pnorm(-q[j*I+0]/sdq0[j], 0, 1, 1, 0))+pnorm(-q[j*I+0]/sdq0[j], 0, 1, 1, 0);
	  propq0=qnorm(temp, 0, 1, 1, 0)*sdq0[j]+q[j*I+0];	
	       
	  // 2.1.2. Calculate the acceptance ratio for q1
	  Cq0=log(pnorm((1-q[j*I+0])/sdq0[j], 0, 1, 1, 0)-pnorm(-q[j*I+0]/sdq0[j],0,1,1,0))-log(pnorm((1-propq0)/sdq0[j], 0, 1, 1, 0)-pnorm(-propq0/sdq0[j],0,1,1,0));// log(PHI((1-q0^t)/sd_q0)-PHI(-q0^t/sd_q0)))-log(PHI((1-q0')/sd_q0)-PHI(-q0'/sd_q0)))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
	  logAq0=n[1][1]*(log(1-propq0*c)-log(1-q[j*I+0]*c))+(n[1][0]+n[0][1]+AQ[j]-1)*(log(propq0)-log(q[j*I+0]))+(n[0][0]+BQ[j]-1)*(log(1-propq0)-log(1-q[j*I+0]))+Cq0;
	 
	  // 2.1.3. sample q0
	  U=runif(0,1);
	  if(log(U)<logAq0)
	    {
	      q[j*I+0]=propq0;
	      rate1[j]=rate1[j]+1.0;
	    }
	  // 2.1.4. claculate q1 by use the q1=1-q0/c
	  q[j*I+1]=1-q[j*I+0]*c;

	  // 2.1.4. claculate logAc1 and logAc2 which is used later when sample the common protion c.  
	  logAc1=logAc1+stateX[j*S+0]*log(1/(1+propc))+(1-stateX[j*S+0])*log(1-1/(1+propc))+n[1][1]*log(1-q[j*I+0]*propc)+n[1][0]*log(q[j*I+1]*propc);
	  logAc2=logAc2+stateX[j*S+0]*log(1/(1+c))+(1-stateX[j*S+0])*log(1-1/(1+c))+n[1][1]*log(1-q[j*I+0]*c)+n[1][0]*log(q[j*I+1]*c);

	  // 2.2 sample inner mixture portion pi
	  a=Api[j]+nZ1[j];
	  b=Bpi[j]+nZ0[j];
	  pi[j]=rbeta(a,b);
	  //printf("j, a, b and pi= %d %lf %lf %lf\n",j, a, b, pi[j]);
	  
	  // 2.3. sample lambda0 and mu1,phi1 when method="ZIP&NB"
	  if (*met==0)
	    {
	      // 2.3.1. sample lambda0
	      a=Alambda[j*I+0]+sumZ1[j];
	      b=1/(Blambda[j*I+0]+nZ1[j]);
	      //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda[0], sumZ1, Blambda[0], nZ1);
	      lambda[j*I+0]=rgamma(a,b);
	      
	      // 2.3.2. sample mu1
	      U1=runif(0,1);
	      temp=U1*pnorm(mu[j*I+1]/sdmu1[j], 0, 1, 1, 0)+pnorm(-mu[j*I+1]/sdmu1[j], 0, 1, 1, 0);
	      propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu1[j]+mu[j*I+1];	  
	      Cmu1=pnorm(mu[j*I+1]/sdmu1[j],0,1,1,1)-pnorm(propmu1/sdmu1[j],0,1,1,1);
	      logAmu1=sum1[j]*(log(propmu1/(propmu1+phi[j*I+1]))-log(mu[j*I+1]/(mu[j*I+1]+phi[j*I+1])))+phi[j*I+1]*n1[j]*(log(mu[j*I+1]+phi[j*I+1])-log(propmu1+phi[j*I+1]))+(Amu[j*I+1]-1.0)*log(propmu1/mu[j*I+1])-Bmu[j*I+1]*(propmu1-mu[j*I+1])+Cmu1;
	      U=runif(0,1);
	      if(log(U)<logAmu1)
		{
		  mu[j*I+1]=propmu1;
		  rate[j*4+0]=rate[j*4+0]+1.0;
		}
	      // 2.3.3. sample phi1
	      U1=runif(0,1);
	      temp=U1*pnorm(phi[j*I+1]/sdphi1[j], 0, 1, 1, 0)+pnorm(-phi[j*I+1]/sdphi1[j], 0, 1, 1, 0);
	      propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi1[j]+phi[j*I+1];	
	      
	      Cphi1=pnorm(phi[j*I+1]/sdphi1[j],0,1,1,1)-pnorm(propphi1/sdphi1[j],0,1,1,1);
	      sumgammaphi1=0;
	      for (s=0;s<S;s++)
		{
		  if ((stateX[j*S+s]==1)&(data[j*S+s]>0))
		    {
		      for (j2=0;j2<data[j*S+s];j2++)
			{
			  sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi[j*I+1]+j2);
			}
		    }
		}
	      logAphi1=sumgammaphi1+sum1[j]*(log(mu[j*I+1]+phi[j*I+1])-log(mu[j*I+1]+propphi1))+n1[j]*(propphi1*log(propphi1/(mu[j*I+1]+propphi1))-phi[j*I+1]*log(phi[j*I+1]/(mu[j*I+1]+phi[j*I+1])))+(Aphi[j*I+1]-1.0)*log(propphi1/phi[j*I+1])-Bphi[j*I+1]*(propphi1-phi[j*I+1])+Cphi1;
	      U=runif(0,1);
	      if(log(U)<logAphi1)
		{
		  phi[j*I+1]=propphi1;
		  rate[j*4+1]=rate[j*4+1]+1.0;
		}
	    }

	  // 2.4. sample lambda when method="poisson"(*met=1)
	  if (*met==1)
	    {
	      a=Alambda[j*I+0]+sumZ1[j];
	      b=1/(Blambda[j*I+0]+nZ1[j]);
	      //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda[0], sumZ1, Blambda[0], nZ1);
	      lambda[j*I+0]=rgamma(a,b);
	      a=Alambda[j*I+1]+sum1[j];
	      b=1/(Blambda[j*I+1]+n1[j]);
	      lambda[j*I+1]=rgamma(a,b);
	      //Rprintf("sample lambda1 %lf %lf %lf %d\n", Alambda[1], sum1, Blambda[1], n1);
	      if(lambda[j*I+0]>lambda[j*I+1])
		{
		  temp=lambda[j*I+1];
		  lambda[j*I+1]=lambda[j*I+0];
		  lambda[j*I+0]=temp;
		  //Rprintf("lambda= %lf %lf\n", lambda[0], lambda[1]);
		}
	      //Rprintf("lambda= %lf %lf\n", lambda[0], lambda[1]);
	    }
	  
	  // 2.5. sample mu and phi when method="NB" (*met=2)
	  if (*met==2) // sample hiarchique parameters phi_0 and phi_1 of Gamma distribution by Metropolis-Hastling method  
	    {
	      // 2.5.1.Sample proposal mu from truncated normal distributions (postive tail)
	      U1=runif(0,1);
	      temp=U1*pnorm(mu[j*I+0]/sdmu0[j], 0, 1, 1, 0)+pnorm(-mu[j*I+0]/sdmu0[j], 0, 1, 1, 0);
	      propmu0=qnorm(temp, 0, 1, 1, 0)*sdmu0[j]+mu[j*I+0];	
	      U1=runif(0,1);
	      temp=U1*pnorm(mu[j*I+1]/sdmu1[j], 0, 1, 1, 0)+pnorm(-mu[j*I+1]/sdmu1[j], 0, 1, 1, 0);
	      propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu1[j]+mu[j*I+1];	
	      
	      // 2.5.2. Calculate the acceptance ratio for mu
	      Cmu0=pnorm(mu[j*I+0]/sdmu0[j],0,1,1,1)-pnorm(propmu0/sdmu0[j],0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
	      Cmu1=pnorm(mu[j*I+1]/sdmu1[j],0,1,1,1)-pnorm(propmu1/sdmu1[j],0,1,1,1);
	      logAmu0=sumZ1[j]*(log(propmu0/(propmu0+phi[j*I+0]))-log(mu[j*I+0]/(mu[j*I+0]+phi[j*I+0])))+phi[j*I+0]*nZ1[j]*(log(mu[j*I+0]+phi[j*I+0])-log(propmu0+phi[j*I+0]))+(Amu[j*I+0]-1.0)*log(propmu0/mu[j*I+0])-Bmu[j*I+0]*(propmu0-mu[j*I+0])+Cmu0;
	      logAmu1=sum1[j]*(log(propmu1/(propmu1+phi[j*I+1]))-log(mu[j*I+1]/(mu[j*I+1]+phi[j*I+1])))+phi[j*I+1]*n1[j]*(log(mu[j*I+1]+phi[j*I+1])-log(propmu1+phi[j*I+1]))+(Amu[j*I+1]-1.0)*log(propmu1/mu[j*I+1])-Bmu[j*I+1]*(propmu1-mu[j*I+1])+Cmu1;
	      // 2.5.3. sample mu
	      U=runif(0,1);
	      if(log(U)<logAmu0)
		{
		  mu[j*I+0]=propmu0;
		  rate[j*4+2]=rate[j*4+2]+1.0;
		}
	      U=runif(0,1);
	      if(log(U)<logAmu1)
		{
		  mu[j*I+1]=propmu1;
		  rate[j*4+0]=rate[j*4+0]+1.0;
		}
	      
	      // 2.5.4. Sample proposal phi from truncated normal distribution (postive tail)
	      U1=runif(0,1);
	      temp=U1*pnorm(phi[j*I+0]/sdphi0[j], 0, 1, 1, 0)+pnorm(-phi[j*I+0]/sdphi0[j], 0, 1, 1, 0);
	      propphi0=qnorm(temp, 0, 1, 1, 0)*sdphi0[j]+phi[j*I+0];	
	      U1=runif(0,1);
	      temp=U1*pnorm(phi[j*I+1]/sdphi1[j], 0, 1, 1, 0)+pnorm(-phi[j*I+1]/sdphi1[j], 0, 1, 1, 0);
	      propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi1[j]+phi[j*I+1];	
	      
	      // 2.5.5. Calculate the acceptance ratio for phi
	      Cphi0=pnorm(phi[j*I+0]/sdphi0[j],0,1,1,1)-pnorm(propphi0/sdphi0[j],0,1,1,1);// log(PHI(phi^t/sd_phi))-log(PHI(phi'/sd_phi))##using R function pnorm(x, mu, phi, iflowtail=1, iflog=1)
	      Cphi1=pnorm(phi[j*I+1]/sdphi1[j],0,1,1,1)-pnorm(propphi1/sdphi1[j],0,1,1,1);
	      sumgammaphi0=0;
	      sumgammaphi1=0;
	      for (s=0;s<S;s++)
		{
		  if ((stateX[j*S+s]==0)&(stateZ[j*S+s]==1)&(data[j*S+s]>0))
		    {
		      for (j2=0;j2<data[j*S+s];j2++)
			{
			  sumgammaphi0=sumgammaphi0+log(propphi0+j2)-log(phi[j*I+0]+j2);
			}
		    }
		  if ((stateX[j*S+s]==1)&(data[j*S+s]>0))
		    {
		      for (j2=0;j2<data[j*S+s];j2++)
			{
			  sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi[j*I+1]+j2);
			}
		    }	
		}
	      logAphi0=sumgammaphi0+sumZ1[j]*(log(mu[j*I+0]+phi[j*I+0])-log(mu[j*I+0]+propphi0))+nZ1[j]*(propphi0*log(propphi0/(mu[j*I+0]+propphi0))-phi[j*I+0]*log(phi[j*I+0]/(mu[j*I+0]+phi[j*I+0])))+(Aphi[j*I+0]-1.0)*log(propphi0/phi[j*I+0])-Bphi[j*I+0]*(propphi0-phi[j*I+0])+Cphi0;
	      logAphi1=sumgammaphi1+sum1[j]*(log(mu[j*I+1]+phi[j*I+1])-log(mu[j*I+1]+propphi1))+n1[j]*(propphi1*log(propphi1/(mu[j*I+1]+propphi1))-phi[j*I+1]*log(phi[j*I+1]/(mu[j*I+1]+phi[j*I+1])))+(Aphi[j*I+1]-1.0)*log(propphi1/phi[j*I+1])-Bphi[j*I+1]*(propphi1-phi[j*I+1])+Cphi1;
	      // 2.5.6. sample phi
	      U=runif(0,1);
	      if(log(U)<logAphi0)
		{
		  phi[j*I+0]=propphi0;
		  rate[j*4+3]=rate[j*4+3]+1.0;
		}
	      U=runif(0,1);
	      if(log(U)<logAphi1)
		{
		  phi[j*I+1]=propphi1;
		  rate[j*4+1]=rate[j*4+1]+1.0;
		}
	      //Rprintf("mu,phi= %lf %lf %lf %lf\n",mu[j*I+1], phi[j*I+1], mu[j*I+0],  phi[j*I+0]);
	      if(mu[j*I+0]>mu[j*I+1])
		{
		  //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu[0], mu[1], phi[0], phi[1]);
		  temp=phi[j*I+1];
		  phi[j*I+1]=phi[j*I+0];
		  phi[j*I+0]=temp;
		  temp=mu[j*I+1];
		  mu[j*I+1]=mu[j*I+0];
		  mu[j*I+0]=temp;
		}
	      //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu[0], mu[1], phi[0], phi[1]);
	    }
	}
            
      // 2.1.0.2. Calculate the acceptance ratio for c
      Cc=pnorm(c/sdc,0,1,1,1)-pnorm(propc/sdc,0,1,1,1);//formula is log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu)), where PHI is standard norm cdf##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
      logAc=logAc1-logAc2+Cc;
      
      // 2.1.0.3. sample c
      U=runif(0,1);
      if(log(U)<logAc)
	{
	  c=propc;
	  rate1[nr]=rate1[nr]+1.0;
	}
	 
      // iteration step 3. calculate likelihood function 
      logl=0.0;  
      for (j=0;j<nr;j++)
	{     
	  pX1[j]=1/(1+exp(log(1-q[j*I+1])-log(q[j*I+0])));//stationary initial distribution
	  pX0[j]=1-pX1[j];
	}
      for (s=0; s<S; s++)
	{
	  for (j=0;j<nr;j++)
	    {     
	      logl=logl+log(((1-pi[j])*indexY[j*S+s]+pi[j]*pY[0][j*S+s]*indexY[j*S+s]+pY[0][j*S+s]*(1-indexY[j*S+s]))*pX0[j]+pY[1][j*S+s]*pX1[j]);
	    }
	}
      //Rprintf("loglikelihood= %lf\n", logl);
      
      // iteration step 4. pass the results back to R
      if ((z>=bN)&(z==j1*10+bN))
	{
	  for (j=0;j<nr;j++)
	    {
	      *(es_q0+j1*nr+j)=q[j*I+0];
	      *(es_q1+j1*nr+j)=q[j*I+1];
	      *(es_pi+j1*nr+j)=pi[j];
	      *(es_lambda1+j1*nr+j)=lambda[j*I+1];
	      *(es_lambda0+j1*nr+j)=lambda[j*I+0];
	      *(es_mu1+j1*nr+j)=mu[j*I+1];
	      *(es_phi1+j1*nr+j)=phi[j*I+1];
	      *(es_mu0+j1*nr+j)=mu[j*I+0];
	      *(es_phi0+j1*nr+j)=phi[j*I+0];
	    }
	  *(loglikeli+j1)=logl; 
	  j1++;
	}
      if (z==z1*jN)
	{
	  Rprintf("%d % ",z1*10);
	  z1++;	    
	}

/*      if (z==j1*1000)
      	{
	  j1++;
	  Rprintf("MCMC step= %d \n", z+1);
	  for (j=0;j<nr;j++)
	    {
	      Rprintf("parameters of %d", j+1);
	      Rprintf("th experiment are %lf %lf %lf %lf %lf\n", mu[j*I+1], phi[j*I+1], pi[j], mu[j*I+0], phi[j*I+0]);
	      Rprintf("and the transition probability and ratio= %lf %lf %lf %lf\n", q[j*I+1], q[j*I+0], c, q[j*I+0]/(q[j*I+0]+1-q[j*I+1]));
	    }
	}*/
    }
  Rprintf("\n");
  PutRNGstate();
  double iter;
  iter=*N;
  
  for (i=0;i<nr; i++)
    {
      *(acrate+i*4+0)=rate[i*4+0]/iter;
      *(acrate+i*4+1)=rate[i*4+1]/iter;
      *(acrate+i*4+2)=rate[i*4+2]/iter;
      *(acrate+i*4+3)=rate[i*4+3]/iter;
      //Rprintf("The acceptance rate for %d", i+1);
      //Rprintf("th experiment are:%lf %lf %lf %lf \n",rate[i*4+0]/iter,rate[i*4+1]/iter, rate[i*4+2]/iter, rate[i*4+3]/iter);// acceptance rate
    }

  //Rprintf("The acceptance rate for transition probabilities are:");
  for (i=0;i<nr;i++)
    { 
      *(acrate1+i)=rate1[i]/iter;
      //Rprintf(" %lf ",rate1[i]/iter);// acceptance rate
    }
 // Rprintf("%lf \n", rate1[nr]/iter);
  *(acrate1+nr)=rate1[nr]/iter;

  for (i=0;i<nr;i++)
    {
      for (s=0; s<S; s++)
	{
	  *(PP+i*S+s)=sumprobX[i*S+s]/iter;
	}
    }
  //Rprintf("iter= %lf\n", iter);
  //for(s=47;s<51;s++)
  //{
  //  Rprintf("sumprobX[47:51,]= %lf %lf\n", sumprobX[s], sumprobX[S+s]);
  //  Rprintf("PP[47:51,]= %lf %lf \n", *(PP+s), *(PP+S+s));
  //}
  
  for (i=0; i<I; i++)
    {
      free(pY[i]);
      free(pXcY[i]);
      free(logpY[i]);
    }
  
  free(pY);
  free(pXcY);
  free(logpY);
  free(probX);
}


