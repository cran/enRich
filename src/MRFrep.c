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
void MRFrep(int *data1, int *nrep, int *size, int *met, double *qprior, double *piprior, double *poisprior, double *NBprior,  int *N, int *Nb, int *jumpN, double *var, double *PP, double *es_pi, double *es_q1, double *es_q0, double *es_lambda1, double *es_mu1, double *es_phi1, double *es_lambda0, double *es_mu0, double *es_phi0, double *loglikeli, double *acrate)
{
  int S=*size;
  int N1=*N;
  int bN=*Nb;
  int jN=*jumpN;
  int nr=*nrep;
  int nrI=nr*I;
  int nrS=nr*S;
  //Rprintf("S, N1, nr, nrI and nrS= %d %d %d %d %d\n", S, N1, nr, nrI, nrS);
  double AQ[I], BQ[I];  // priors of transition probabilities.
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
  stateX=(int*)malloc(sizeof(int*)*S);
  int *stateZ; // inner latent variables
  stateZ=(int*)malloc(sizeof(int*)*nrS);
  double *lambda, *mu, *phi, q[I], *pi; // parameters of Poisson-Gamma mixing distribution and transition probabilities
  lambda=(double*)malloc(sizeof(double*)*nrI);
  mu=(double*)malloc(sizeof(double*)*nrI);
  phi=(double*)malloc(sizeof(double*)*nrI);
  pi=(double*)malloc(sizeof(double*)*nr);
  int i, j, j1, j2, k, l, s, z; 
 
  double *probX, *sumprobX, *probZ, pX1[I];
  int *indexY;
  probX=(double*)malloc(sizeof(double*)*S);
  sumprobX=(double*)malloc(sizeof(double*)*S);
  probZ=(double*)malloc(sizeof(double*)*nr);
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
      pXcY[i]=(double*)malloc(sizeof(double*)*S);
    }
  
  double *sum2X;//using for counting pairs of X
  sum2X=(double*)malloc(sizeof(double*)*(S-1));
  double U,U1,temp,a,b,logl;
  
  ///////// R data1 pass to data
  for(s=0;s<nrS;s++)
    {
      data[s]=*(data1+s);
    }
  for (s=0;s<S;s++)
    {
      sumprobX[s]=0.0;
    }
  //for(s=48605;s<48610;s++)
  //{
  //  Rprintf("data[48606:48610,]= %d %d\n", data[s], data[S+s]);
  //}

  ///////// give priors for parameters; 
  AQ[0]=*(qprior+0);
  AQ[1]=*(qprior+1);
  BQ[0]=*(qprior+2);
  BQ[1]=*(qprior+3);  
  //Rprintf("prior of q= %f %f %f %f \n",AQ[0], AQ[1], BQ[0], BQ[1]);
  
  for (i=0;i<nr;i++)
    {
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
      //Rprintf("Prior of %d",i);
      //Rprintf("th experiment, prior of pi= %f %f \n", Api[i], Bpi[i]);
      //Rprintf("prior of lambda= %f %f %f %f \n", Alambda[i*I+1], Blambda[i*I+1], Alambda[i*I+0], Blambda[i*I+0]);
      //Rprintf("prior of mu and phi= %f %f %f %f %f %f %f %f \n", Amu[i*I+1], Bmu[i*I+1], Aphi[i*I+1], Bphi[i*I+1], Amu[i*I+0], Bmu[i*I+0], Aphi[i*I+0], Bphi[i*I+0]); 
    }
  
  //////// give initial values, set X=1 if Y>cp, where cp is initial cutting points
  double *sumZ1, *sum1;
  sumZ1=(double*)malloc(sizeof(double*)*nr);
  sum1=(double*)malloc(sizeof(double*)*nr);
  int n1,n0, *nZ1,*nZ0;
  nZ1=(int*)malloc(sizeof(int*)*nr);
  nZ0=(int*)malloc(sizeof(int*)*nr);

  n1=0;
  n0=0;
  for (s=0;s<S;s++)
    {
      temp=0;
      for (i=0;i<nr;i++)
	{
	  temp=temp+data[i*S+s];
	}
      if (temp>cp*(nr-1))
	{
	  stateX[s]=1;
	  n1=n1++;
	}
      else
	{
	  stateX[s]=0;
	  n0=n0++;
	}
    }
  
  for (i=0;i<nr;i++)
    {
      sumZ1[i]=0.0;
      sum1[i]=0.0;
      nZ1[i]=0;
      nZ0[i]=0;
      for (s=0;s<S;s++)
	{ 
	  if(stateX[s]==1)
	    { 
	      sum1[i]=sum1[i]+data[i*S+s];
	      stateZ[i*S+s]=2;
	      indexY[i*S+s]=0;
	    }
	  else
	    {
	      if (data[i*S+s]==0)
		{
		  indexY[i*S+s]=1;
		  stateZ[i*S+s]=0;
		  nZ0[i]=nZ0[i]++;
		}
	      else
		{
		  indexY[i*S+s]=0;
		  stateZ[i*S+s]=1;
		  sumZ1[i]=sumZ1[i]+data[i*S+s];
		  nZ1[i]=nZ1[i]++;
		}
	    }
	}
      if (sum1[i]==0|n1==0)
	{
	  sum1[i]=5;
	  n1=1;
	}
      //printf("initial sumZ1, sum1, sum, n1, n0, n1+n0=% lf %lf %lf %d %d %d \n",sumZ1, sum1, sum, nZ0, nZ1, nZ0+nZ1, n0, n1, n0+n1);
      if (*met==0)
	{
	  lambda[i*I+0]=sumZ1[i]/nZ1[i];
	  lambda[i*I+1]=0;
	  mu[i*I+1]=sum1[i]/n1;
	  mu[i*I+0]=0;
	  phi[i*I+1]=1;
	  phi[i*I+0]=0;
	}
      if (*met==1)
	{
	  lambda[i*I+0]=sumZ1[i]/nZ1[i];//sum0/n0;
	  lambda[i*I+1]=sum1[i]/n1; 
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
	  mu[i*I+1]=sum1[i]/n1;
	  phi[i*I+0]=1;
	  phi[i*I+1]=1;
	  //printf("initial mu and pi %lf %lf %lf\n", mu[0],mu[1], phi[0], phi[1]);
	}
      pi[i]=(nZ1[i]+0.0)/(nZ0[i]+0.0);
      //printf("initial nZ and pi %d %d %lf\n", nZ1, nZ0, pi);
    }
  
  //for(s=48605;s<48610;s++)
  //{
  //  Rprintf("indexY[48606:48610,]= %d %d\n", indexY[s], indexY[S+s]);
  //}
  
  //Rprintf("initial values\n");
  //Rprintf("transition probability= %lf %lf\n", q[1], q[0]);
  //for (j=0;j<nr;j++)
  //{
  //  Rprintf("parameters of %d", j+1);
  //  Rprintf("th experiment are %lf %lf %lf %lf %lf\n", pi[j], mu[j*I+1], phi[j*I+1], mu[j*I+0], phi[j*I+0]);
  //}
  
  q[0]=0.1; //prob(X_s=1|X_{s-1}=0)
  q[1]=0.9; //prob(X_s=1|X_{s-1}=1)

  double *sdmu0, *sdmu1, *sdphi0, *sdphi1, *rate, sumgammaphi0, sumgammaphi1;
  sdmu0=(double*)malloc(sizeof(double*)*nr);
  sdmu1=(double*)malloc(sizeof(double*)*nr);
  sdphi0=(double*)malloc(sizeof(double*)*nr); 
  sdphi1=(double*)malloc(sizeof(double*)*nr); 
  int nr4=4*nr;
  rate=(double*)malloc(sizeof(double*)*nr4);
  for (i=0;i<nr;i++)
    {
      sdmu1[i]=*(var+i*4+0);
      sdphi1[i]=*(var+i*4+1);
      sdmu0[i]=*(var+i*4+2);
      sdphi0[i]=*(var+i*4+3);
      rate[i*4+0]=0.0;
      rate[i*4+1]=0.0;
      rate[i*4+2]=0.0;
      rate[i*4+3]=0.0;
      //Rprintf("var= %lf %lf %lf %lf \n", sdmu1[i], sdphi1[i], sdmu0[i], sdphi0[i]);
    }
  double propmu0, propmu1, propphi0, propphi1, Cmu0, Cmu1, logAmu0, logAmu1, Cphi0, Cphi1, logAphi0, logAphi1;
 
  j1=0;
  int z1=0;
  GetRNGstate();
  for(z=0;z<*N;z++)//show results for every 10th iterations
    {      
      Q[0][0]=1-q[0];
      Q[0][1]=q[0];
      Q[1][0]=1-q[1];
      Q[1][1]=q[1];
      n[0][0]=0;
      n[0][1]=0;
      n[1][0]=0;
      n[1][1]=0;
      n0=0;
      n1=0;

      // iteration step 1: Given parameters renew the states X of bins by using direct gibbs method. We also sample the inner states Z given X==0. 
      // 1.1 sample stateX
      for (s=0; s<S; s++)
	{
	  for (j=0;j<nr;j++)
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
	}
      //Rprintf("mu and phi are: %lf %lf %lf %lf %lf %lf %lf %lf\n", mu[1], phi[1], mu[0], phi[0], mu[3], phi[3], mu[2], mu[2]);
      //for(s=48605;s<48610;s++)
      //{
      //  Rprintf("pY[48606:48610,]= %lf %lf %lf %lf \n", pY[1][s], pY[0][s], pY[1][S+s], pY[0][S+s]);
      //}
      pXcY[0][0]=1.0;
      pXcY[1][0]=1.0;
      for (j=0;j<nr;j++)
	{
	  pXcY[0][0]=pXcY[0][0]*((1-pi[j])*indexY[j*S+0]+pi[j]*indexY[j*S+0]*pY[0][j*S+0]+(1-indexY[j*S+0])*pY[0][j*S+0]);
	  pXcY[1][0]=pXcY[1][0]*pY[1][j*S+0];
	}
      pXcY[0][0]=pXcY[0][0]*Q[0][stateX[1]];
      pXcY[1][0]=pXcY[1][0]*Q[1][stateX[1]];
      probX[0]=1/(1+exp(log(pXcY[0][0])-log(pXcY[1][0])));//posterior distribution of P(X_1=1)
      sumprobX[0]=sumprobX[0]+probX[0];
      U=runif(0,1);      
      if(U<probX[0])
	{
	  stateX[0]=1;
	  n1=n1++;//To get sum_s I(X_s=1)
	}
      else
	{
	  stateX[0]=0;
	  n0=n0++;
	}
      for (s=1; s<S-1; s++)
	{
	  pXcY[0][s]=1.0;
	  pXcY[1][s]=1.0;
	  for (j=0;j<nr;j++)
	    {
	      pXcY[0][s]=pXcY[0][s]*((1-pi[j])*indexY[j*S+s]+pi[j]*indexY[j*S+s]*pY[0][j*S+s]+(1-indexY[j*S+s])*pY[0][j*S+s]);
	      pXcY[1][s]=pXcY[1][s]*pY[1][j*S+s];
	    }
	  pXcY[0][s]=pXcY[0][s]*Q[stateX[s-1]][0]*Q[0][stateX[s+1]];
	  pXcY[1][s]=pXcY[1][s]*Q[stateX[s-1]][1]*Q[1][stateX[s+1]];
	  probX[s]=1/(1+exp(log(pXcY[0][s])-log(pXcY[1][s])));//posterior distribution of P(X_1=1)
	  sumprobX[s]=sumprobX[s]+probX[s];
	  U=runif(0,1);
	  if(U<probX[s])
	    {
	      stateX[s]=1;
	      n1=n1++;//To get sum_s I(X_s=1)
	    }     
	  else
	    {
	      stateX[s]=0; 
	      n0=n0++;
	    }
	     
	}
      pXcY[0][S-1]=1.0;
      pXcY[1][S-1]=1.0;
      for (j=0;j<nr;j++)
	{
	  pXcY[0][S-1]=pXcY[0][S-1]*((1-pi[j])*indexY[j*S+S-1]+pi[j]*indexY[j*S+S-1]*pY[0][j*S+S-1]+(1-indexY[j*S+S-1])*pY[0][j*S+S-1]);
	  pXcY[1][S-1]=pXcY[1][S-1]*pY[1][j*S+S-1];
	}
      pXcY[0][S-1]=pXcY[0][S-1]*Q[stateX[S-2]][0];
      pXcY[1][S-1]=pXcY[1][S-1]*Q[stateX[S-2]][1];
      probX[S-1]=1/(1+exp(log(pXcY[0][S-1])-log(pXcY[1][S-1])));//posterior distribution of P(X_1=1)
      sumprobX[S-1]=sumprobX[S-1]+probX[S-1];
      U=runif(0,1);
      if(U<probX[S-1])
	{
	  stateX[S-1]=1;
	  n1=n1++;//To get sum_s I(X_s=1)
	}     
      else
	{	  
	  stateX[S-1]=0;
	  n0=n0++;
	}
      
      // 1.2 given stateX we sample stateZ
      for (j=0;j<nr;j++)
	{
	  sumZ1[j]=0.0;
	  sum1[j]=0.0;
	  nZ0[j]=0;
	  nZ1[j]=0;
	  if (*met<2)
	    {
	      probZ[j]=1/(1+exp(log(1-pi[j])-log(pi[j])+lambda[j*I+0]));
	    }
	  if (*met==2)
	    {
	      probZ[j]=1/(1+exp(log(1-pi[j])-log(pi[j])-phi[j*I+0]*log(phi[j*I+0]/(phi[j*I+0]+mu[j*I+0]))));
	    }
	  //Rprintf("probZ= %lf\n", probZ);
	  for (s=0;s<S;s++)
	    {
	      if (stateX[s]==0)
		{
		  if (data[j*S+s]>0)
		    {
		      stateZ[j*S+s]=1;
		      sumZ1[j]=sumZ1[j]+data[j*S+s];
		      nZ1[j]=nZ1[j]++;
		    }
		  else
		    {
		      U=runif(0,1);
		      if (U<probZ[j])
			{
			  stateZ[j*S+s]=1;
			  sumZ1[j]=sumZ1[j]+data[j*S+s]; //To get sum_s Y_sI(X_s=0, Z_s=1)
			  nZ1[j]=nZ1[j]++;//To get sum_s I(X_s=0, Z_s=1)
			}
		      else
			{
			  stateZ[j*S+s]=0;
			  nZ0[j]=nZ0[j]++;//To get sum_s I(X_s=0, Z_s=0)
			}
		    }
		}
	      else
		{		  
		  sum1[j]=sum1[j]+data[j*S+s]; //To get sum_s Y_sI(X_s=1)
		  stateZ[j*S+s]=2;
		}
	    }
	}
    
      // iteration step 2: sample the parameters from posterior distribution.
      // 2.1 sample transition probabilities from posterior beta distribution
      // 2.1.1 count the pair
      for (s=0;s<S-1;s++)
	{
	  sum2X[s]=stateX[s]+stateX[s+1]*2;
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
      
      // 2.1.2 sample transition probabilities
      
      for(i=0;i<I;i++)
	{ 
	  a=AQ[i]+n[i][1];
	  b=BQ[i]+n[i][0];
	  q[i]=rbeta(a,b);
	  //Rprintf("sample q %lf %lf %lf\n", a, b, q[i]);
	}
	
      // 2.2 sample inner mixture portion pi
      for (j=0;j<nr;j++)
	{
	  a=Api[j]+nZ1[j];
	  b=Bpi[j]+nZ0[j];
	  pi[j]=rbeta(a,b);
	  //printf("pi= %lf\n", pi[j]);
	  //pi=1;//no-zero inflated model
	  

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
	      logAmu1=sum1[j]*(log(propmu1/(propmu1+phi[j*I+1]))-log(mu[j*I+1]/(mu[j*I+1]+phi[j*I+1])))+phi[j*I+1]*n1*(log(mu[j*I+1]+phi[j*I+1])-log(propmu1+phi[j*I+1]))+(Amu[j*I+1]-1.0)*log(propmu1/mu[j*I+1])-Bmu[j*I+1]*(propmu1-mu[j*I+1])+Cmu1;
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
		  if (stateX[s]==1&data[j*S+s]>0)
		    {
		      for (j2=0;j2<data[j*S+s];j2++)
			{
			  sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi[j*I+1]+j2);
			}
		    }
		}
	      logAphi1=sumgammaphi1+sum1[j]*(log(mu[j*I+1]+phi[j*I+1])-log(mu[j*I+1]+propphi1))+n1*(propphi1*log(propphi1/(mu[j*I+1]+propphi1))-phi[j*I+1]*log(phi[j*I+1]/(mu[j*I+1]+phi[j*I+1])))+(Aphi[j*I+1]-1.0)*log(propphi1/phi[j*I+1])-Bphi[j*I+1]*(propphi1-phi[j*I+1])+Cphi1;
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
	      b=1/(Blambda[j*I+1]+n1);
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
	      logAmu1=sum1[j]*(log(propmu1/(propmu1+phi[j*I+1]))-log(mu[j*I+1]/(mu[j*I+1]+phi[j*I+1])))+phi[j*I+1]*n1*(log(mu[j*I+1]+phi[j*I+1])-log(propmu1+phi[j*I+1]))+(Amu[j*I+1]-1.0)*log(propmu1/mu[j*I+1])-Bmu[j*I+1]*(propmu1-mu[j*I+1])+Cmu1;
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
		  if (stateX[s]==0&stateZ[j*S+s]==1&data[j*S+s]>0)
		    {
		      for (j2=0;j2<data[j*S+s];j2++)
			{
			  sumgammaphi0=sumgammaphi0+log(propphi0+j2)-log(phi[j*I+0]+j2);
			}
		    }
		  if (stateX[s]==1&data[j*S+s]>0)
		    {
		      for (j2=0;j2<data[j*S+s];j2++)
			{
			  sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi[j*I+1]+j2);
			}
		    }	
		}
	      logAphi0=sumgammaphi0+sumZ1[j]*(log(mu[j*I+0]+phi[j*I+0])-log(mu[j*I+0]+propphi0))+nZ1[j]*(propphi0*log(propphi0/(mu[j*I+0]+propphi0))-phi[j*I+0]*log(phi[j*I+0]/(mu[j*I+0]+phi[j*I+0])))+(Aphi[j*I+0]-1.0)*log(propphi0/phi[j*I+0])-Bphi[j*I+0]*(propphi0-phi[j*I+0])+Cphi0;
	      logAphi1=sumgammaphi1+sum1[j]*(log(mu[j*I+1]+phi[j*I+1])-log(mu[j*I+1]+propphi1))+n1*(propphi1*log(propphi1/(mu[j*I+1]+propphi1))-phi[j*I+1]*log(phi[j*I+1]/(mu[j*I+1]+phi[j*I+1])))+(Aphi[j*I+1]-1.0)*log(propphi1/phi[j*I+1])-Bphi[j*I+1]*(propphi1-phi[j*I+1])+Cphi1;
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
      
      // iteration step 3. calculate likelihood function 
      logl=0.0;  
      pX1[1]=1/(1+exp(log(1-q[1])-log(q[0])));//stationary initial distribution
      pX1[0]=1-pX1[1];
      for (s=0; s<S; s++)
	{
	  for (j=0;j<nr;j++)
	    {
	      logl=logl+log(((1-pi[j])*indexY[j*S+s]+pi[j]*pY[0][j*S+s]*indexY[j*S+s]+pY[0][j*S+s]*(1-indexY[j*S+s]))*pX1[0]+pY[1][j*S+s]*pX1[1]);
	    }
	}
      //Rprintf("loglikelihood= %lf\n", logl);
      
      // iteration step 4. pass the results back to R
      //Rprintf("%d %d\n", z, j1*10+bN);
      if (z>=bN&z==j1*10+bN)
	{
	  *(es_q0+j1)=q[0];
	  *(es_q1+j1)=q[1];
	  for (j=0;j<nr;j++)
	    {
	      *(es_pi+j1*nr+j)=pi[j];
	      *(es_lambda1+j1*nr+j)=lambda[j*I+1];
	      *(es_lambda0+j1*nr+j)=lambda[j*I+0];
	      *(es_mu1+j1*nr+j)=mu[j*I+1];
	      *(es_phi1+j1*nr+j)=phi[j*I+1];
	      *(es_mu0+j1*nr+j)=mu[j*I+0];
	      *(es_phi0+j1*nr+j)=phi[j*I+0];
	    }
	  *(loglikeli+j1)=logl;  
	  j1=j1++;
	}    
      if (z==z1*jN)
	{
	   Rprintf("%d % ",z1*10);
	  z1++;	    
	}      
      
      /* if (z==j1*1000)
	 {
	 j1=j1++;
	 Rprintf("MCMC step= %d \n", z);
	 Rprintf("transition probability= %lf %lf\n", q[1], q[0]);
	 for (j=0;j<nr;j++)
	 {
	 Rprintf("parameters of %d", j+1);
	 Rprintf("th experiment are %lf %lf %lf %lf %lf\n", mu[j*I+1], phi[j*I+1], pi[j], mu[j*I+0], phi[j*I+0]);
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
      //Rprintf("The acceptance rate for %d", i);
      //Rprintf("th experiment are:%lf %lf %lf %lf \n",rate[i*4+0]/iter,rate[i*4+1]/iter, rate[i*4+2]/iter, rate[i*4+3]/iter);// acceptance rate
    }

  for (s=0; s<S; s++)
    {
      *(PP+s)=sumprobX[s]/iter;
    }
  
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


