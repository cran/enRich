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
#define I 2 // dimension of transition probabilities


//**********************************************************//
//*********************   Main program  ********************//
void MRF(int *data1, int *size, int *met, double *qprior, double *piprior, double *poisprior, double *NBprior, int *N, int *Nb, int *jumpN, double *var, double *PP, double *es_pi, double *es_q1, double *es_q0, double *es_lambda1, double *es_mu1, double *es_phi1, double *es_lambda0, double *es_mu0, double *es_phi0, double *loglikeli, double *acrate)
{
  int S=*size;
//  int N1=*N;
  int bN=*Nb;
  int jN=*jumpN;
  double AQ[I], BQ[I];  // priors of transition probabilities.
  double Api, Bpi;// priors of mixture portion of ZIP model.
  double Alambda[I], Blambda[I]; // priors of lambda for poisson fitting.
  double Amu[I], Bmu[I];//priors of mu_i, which is the mean of NB for ith component;
  double Aphi[I], Bphi[I];//priors of phi_i, which is the overdispersion of NB for ithe component;
  int *data;// bin counts--observations.
  data=(int*)malloc(sizeof(int*)*S);
  int *stateX; // bin states, latent variables
  stateX=(int*)malloc(sizeof(int*)*S);
  int *stateZ; // inner latent variables
  stateZ=(int*)malloc(sizeof(int*)*S);
  double lambda[I], mu[I],phi[I], q[I], pi; // parameters of Poisson-Gamma mixing distribution and transition probabilities
  phi[0]=0.0;
  phi[1]=0.0;
  mu[0]=0.0;
  mu[1]=0.0;
  lambda[0]=0.0;
  lambda[1]=0.0;
  int i, j1, j2, s, z; 
 
  double *probX, *sumprobX;
  double probZ=0.0;
  double pX1[I];
  int *indexY;
  probX=(double*)malloc(sizeof(double*)*S);
  sumprobX=(double*)malloc(sizeof(double*)*S);
  indexY=(int*)malloc(sizeof(int*)*S);

  int n[I][I];//the counts of pair
  double Q[I][I],**logpY,**pY,**pXcY; 

  logpY=(double**)malloc(sizeof(double*)*I);
  pY=(double**)malloc(sizeof(double*)*I);
  pXcY=(double**)malloc(sizeof(double*)*I);
 
  for(i=0;i<I;i++)
    {
      logpY[i]=(double*)malloc(sizeof(double*)*S);
      pY[i]=(double*)malloc(sizeof(double*)*S);
      pXcY[i]=(double*)malloc(sizeof(double*)*S);
    }
  
  double *sum2X;//using for counting pairs of X
  sum2X=(double*)malloc(sizeof(double*)*(S-1));
  double U, U1,temp,a,b,logl;
  
  ///////// R data1 pass to data
  for(s=0;s<S;s++)
    {
      data[s]=*(data1+s);
      sumprobX[s]=0.0;
    }
  // for(s=971900;s<972000;s++)
  //{
  //  Rprintf("%d ", data[s]);
  //}

  ///////// give priors for parameters; 
  AQ[0]=*(qprior+0);
  AQ[1]=*(qprior+1);
  BQ[0]=*(qprior+2);
  BQ[1]=*(qprior+3);  
  
  Api=*(piprior+0);
  Bpi=*(piprior+1);
  
  Alambda[1]=*(poisprior+0);
  Blambda[1]=*(poisprior+1);
  Alambda[0]=*(poisprior+2);
  Blambda[0]=*(poisprior+3);

  Amu[1]=*(NBprior+0);
  Bmu[1]=*(NBprior+1);
  Aphi[1]=*(NBprior+2);
  Bphi[1]=*(NBprior+3);
  Amu[0]=*(NBprior+4);
  Bmu[0]=*(NBprior+5);
  Aphi[0]=*(NBprior+6);
  Bphi[0]=*(NBprior+7);
  
  //Rprintf("prior of q= %f %f %f %f \n",AQ[0], AQ[1], BQ[0], BQ[1]);
  //Rprintf("prior of pi= %f %f \n", Api, Bpi);
  //Rprintf("prior of lambda= %f %f %f %f \n", Alambda[0], Blambda[0], Alambda[1], Blambda[1]);
  //Rprintf("prior of alpha and beta= %f %f %f %f %f %f %f %f \n", Aalpha[0], Balpha[0], Abeta[0], Bbeta[0], Aalpha[1], Balpha[1], Abeta[1], Bbeta[1]); 
  
  //printf("%d\n", S);
  
  //////// give initial values, set X=1 if Y>cp, where cp is initial cutting points
  double sumZ1, sum1;
  sumZ1=0.0;
  sum1=0.0;
  int n1,n0, nZ1,nZ0;
  n1=0;
  n0=0;
  nZ1=0;
  nZ0=0;
  for (s=0;s<S;s++)
    { 
      if(data[s]>cp)
	{ 
	  stateX[s]=1;
	  sum1=sum1+data[s];
	  n1++;
	  stateZ[s]=2;
	  indexY[s]=0;
	}
      else
	{
	  stateX[s]=0;
	  n0++;
	  if (data[s]==0)
	    {
	      indexY[s]=1;
	      stateZ[s]=0;
	      nZ0++;
	    }
	  else
	    {
	      indexY[s]=0;
	      stateZ[s]=1;
	      sumZ1=sumZ1+data[s];
	      nZ1++;
	    }
	}
    }
  if ((sum1==0)|(n1==0))
    {
      sum1=5;
      n1=1;
    }
  //printf("initial sumZ1, sum1, sum, n1, n0, n1+n0=% lf %lf %lf %d %d %d \n",sumZ1, sum1, sum, nZ0, nZ1, nZ0+nZ1, n0, n1, n0+n1);
  if (*met==0)
    {
      lambda[0]=sumZ1/nZ1;
      lambda[1]=0;
      mu[1]=sum1/n1;
      mu[0]=0;
      phi[1]=1;
      phi[0]=0;
    }
  if (*met==1)
    {
      lambda[0]=sumZ1/nZ1;//sum0/n0;
      lambda[1]=sum1/n1; 
      mu[0]=0;
      mu[1]=0;
      phi[0]=0;
      phi[1]=0;
      //printf("initial lambda %lf %lf\n", lambda[0],lambda[1]); 
    }
  if (*met==2)
    {
      lambda[0]=0;//sum0/n0;
      lambda[1]=0; 
      mu[0]=sumZ1/nZ1;
      mu[1]=sum1/n1;
      phi[0]=1;
      phi[1]=1;
      //printf("initial mu and pi %lf %lf %lf\n", mu[0],mu[1], phi[0], phi[1]);
    }
  pi=(nZ1+0.0)/(nZ0+0.0);
  //printf("initial nZ and pi %d %d %lf\n", nZ1, nZ0, pi);
  
  //j1=48610;
  //Rprintf("Initially, we have data[48610:48630]= %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", data[j1], data[j1+1], data[j1+2], data[j1+3], data[j1+4], data[j1+5], data[j1+6], data[j1+7], data[j1+8], data[j1+9], data[j1+10], data[j1+11], data[j1+12], data[j1+13], data[j1+14], data[j1+15], data[j1+16], data[j1+17], data[j1+18], data[j1+19]);
  //Rprintf("and relative probX= %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", probX[j1], probX[j1+1], probX[j1+2], probX[j1+3], probX[j1+4], probX[j1+5], probX[j1+6], probX[j1+7], probX[j1+8], probX[j1+9], probX[j1+10], probX[j1+11], probX[j1+12], probX[j1+13], probX[j1+14], probX[j1+15], probX[j1+16], probX[j1+17], probX[j1+18], probX[j1+19]);
  //Rprintf("and relative stateX= %d %d %d %d %d  %d %d %d %d %d %d %d %d %d %d  %d %d %d %d %d\n", stateX[j1], stateX[j1+1], stateX[j1+2], stateX[j1+3], stateX[j1+4], stateX[j1+5], stateX[j1+6], stateX[j1+7], stateX[j1+8], stateX[j1+9], stateX[j1+10], stateX[j1+11], stateX[j1+12], stateX[j1+13], stateX[j1+14], stateX[j1+15], stateX[j1+16], stateX[j1+17], stateX[j1+18], stateX[j1+19]);


  //temp=0.0;
  //for (s=0;s<S;s++)
  //  {
  //    temp=temp+indexY[s];
  //  }
  //printf("sum indexY %lf\n", temp);

  q[0]=0.1; //prob(X_s=1|X_{s-1}=0)
  q[1]=0.9; //prob(X_s=1|X_{s-1}=1)

  double sdmu0, sdmu1, sdphi0, sdphi1, sumgammaphi0, sumgammaphi1;
  double propmu0, propmu1, propphi0, propphi1, Cmu0, Cmu1, logAmu0, logAmu1, Cphi0, Cphi1, logAphi0, logAphi1;
  
  sdmu1=*(var+0);
  sdphi1=*(var+1);
  sdmu0=*(var+2);
  sdphi0=*(var+3);
  //Rprintf("var= %lf %lf %lf %lf \n", sdmu1,  sdphi1, sdmu0, sdphi0);
  
  double rate[4];
  rate[0]=0.0;
  rate[1]=0.0;
  rate[2]=0.0;
  rate[3]=0.0; 

  j1=0;
  int z1=0;
  //for(z=0;z<100;z++)
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
      sumZ1=0.0;
      sum1=0.0;
      n0=0;
      n1=0;
      nZ0=0;
      nZ1=0;

      // iteration step 1: Given parameters renew the states X of bins by using direct gibbs method. We also sample the inner states Z given X==0. 
      // 1.1. sample stateX
      for (s=0; s<S; s++)
	{
	  if (*met==0)// poissonNB
	    {
	      logpY[0][s]=dpois(data[s], lambda[0], 1);//use R function dpois(y, lambda, logistrue=1);
	      pY[0][s]=dpois(data[s], lambda[0], 0);//use R function dpois(y, lambda, logisfalse=0);
	      logpY[1][s]=dnbinom(data[s], phi[1], phi[1]/(mu[1]+phi[1]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
	      pY[1][s]=dnbinom(data[s], phi[1], phi[1]/(mu[1]+phi[1]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
	    }
	  if (*met==1)// poisson
	    {
	      for (i=0;i<I;i++)
		{
		  logpY[i][s]=dpois(data[s], lambda[i], 1);//use R function dpois(y, lambda, logistrue=1);
		  pY[i][s]=dpois(data[s], lambda[i], 0);//use R function dpois(y, lambda, logisfalse=0);
		}
	    }
	  if (*met==2)// NB
	    {
	      for (i=0;i<I;i++)
		{
		  logpY[i][s]=dnbinom(data[s], phi[i], phi[i]/(mu[i]+phi[i]), 1);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logistrue=1)
		  pY[i][s]=dnbinom(data[s], phi[i], phi[i]/(mu[i]+phi[i]), 0);//use R function dnbinom(y, size=phi, p=phi/(mu+phi), logisfalse=0)
		}
	    }
	  //Rprintf("%lf ", pY[0][s]*indexY[s]);
	}
      pXcY[0][0]=((1-pi)*indexY[0]+pi*indexY[0]*pY[0][0]+(1-indexY[0])*pY[0][0])*Q[0][stateX[1]];
      pXcY[1][0]=pY[1][0]*Q[1][stateX[1]];
	 if ((pXcY[0][0]<=0)&&(pXcY[1][0]<=0)&&(data[0]>cp))
	{ 
		probX[0]=1;
	}
	else if (pXcY[1][0]<=0)
	{
		probX[0]=0;
	}
	else if (pXcY[0][0]<=0)
	{
		probX[0]=1;
	}
	else
	{
      		probX[0]=1/(1+exp(log(pXcY[0][0])-log(pXcY[1][0])));//posterior distribution of P(X_1=1)
        }
	 sumprobX[0]=sumprobX[0]+probX[0];
      U=runif(0,1);      
      if(U<probX[0])
	{
	  stateX[0]=1;
	  sum1=sum1+data[0]; //To get sum_s Y_sI(X_s=1)
	  n1++;//To get sum_s I(X_s=1)
	}
      else
	{
	  stateX[0]=0;
	  n0++;
	}
      for (s=1; s<S-1; s++)
	{
	  pXcY[0][s]=((1-pi)*indexY[s]+pi*indexY[s]*pY[0][s]+(1-indexY[s])*pY[0][s])*Q[stateX[s-1]][0]*Q[0][stateX[s+1]];
	  pXcY[1][s]=pY[1][s]*Q[stateX[s-1]][1]*Q[1][stateX[s+1]];
	  if ((pXcY[0][s]<=0)&&(pXcY[1][s]<=0)&&(data[s]>cp))
	  { 
		probX[s]=1;
	  }
	  else if (pXcY[1][s]<=0)
	  {
		probX[s]=0;
	  }
	  else if (pXcY[0][s]<=0)
	  {
		probX[s]=1;
	  }
	  else
	  {
      		probX[s]=1/(1+exp(log(pXcY[0][s])-log(pXcY[1][s])));//posterior distribution of P(X_1=1)
	  }
	  sumprobX[s]=sumprobX[s]+probX[s];
	  U=runif(0,1);
	  if(U<probX[s])
	    {
	      stateX[s]=1;
	      sum1=sum1+data[s]; //To get sum_s Y_sI(X_s=1)
	      n1++;//To get sum_s I(X_s=1)
	    }     
	  else
	    {
	      stateX[s]=0; 
	      n0++;
	    }
	     
	}
      pXcY[0][S-1]=((1-pi)*indexY[S-1]+pi*indexY[S-1]*pY[0][S-1]+(1-indexY[S-1])*pY[0][S-1])*Q[stateX[S-2]][0];
      pXcY[1][S-1]=pY[1][S-1]*Q[stateX[S-2]][1];
      if ((pXcY[0][S-1]<=0)&&(pXcY[1][S-1]<=0)&&(data[S-1]>cp))
	{ 
		probX[S-1]=1;
	}
	else if (pXcY[1][S-1]<=0)
	{
		probX[S-1]=0;
	}
	else if (pXcY[0][S-1]<=0)
	{
		probX[S-1]=1;
	}
	else
	{
      	probX[S-1]=1/(1+exp(log(pXcY[0][S-1])-log(pXcY[1][S-1])));//posterior distribution of P(X_1=1)
	}
      sumprobX[S-1]=sumprobX[S-1]+probX[S-1];
      U=runif(0,1);
      if(U<probX[S-1])
	{
	  stateX[S-1]=1;
	  sum1=sum1+data[S-1]; //To get sum_s Y_sI(X_s=1)
	  n1++;//To get sum_s I(X_s=1)
	}     
      else
	{	  
	  stateX[S-1]=0;
	  n0++;
	}
      // 1.2 given stateX sample stateZ
      if (*met<2)
	{
	  probZ=1/(1+exp(log(1-pi)-log(pi)+lambda[0]));
	}
      if (*met==2)
	{
	  probZ=1/(1+exp(log(1-pi)-log(pi)-phi[0]*log(phi[0]/(phi[0]+mu[0]))));
	}
      //Rprintf("probZ= %lf\n", probZ);
      for (s=0;s<S;s++)
	{
	  if (stateX[s]==0)
	    {
	      if (data[s]>0)
		{
		  stateZ[s]=1;
		  sumZ1=sumZ1+data[s];
		  nZ1++;
		}
	      else
		{
		  U=runif(0,1);
		  if (U<probZ)
		    {
		      stateZ[s]=1;
		      sumZ1=sumZ1+data[s]; //To get sum_s Y_sI(X_s=0, Z_s=1)
		      nZ1++;//To get sum_s I(X_s=0, Z_s=1)
		    }
		  else
		    {
		      stateZ[s]=0;
		      nZ0++;//To get sum_s I(X_s=0, Z_s=0)
		    }
		}
	    }
	  else
	    {
	      stateZ[s]=2;
	    }
	}
      //Rprintf("After %d", z);
      //Rprintf("th sampling X and Z,\n");
      //Rprintf("we have sumZ1, sum1 and sumZ1+sum1= %lf %lf %lf \n",sumZ1, sum1, sumZ1+sum1);
      //Rprintf("and nZ0, nZ1, nZ1+nZ0, n0, n1, n0+n1= %d %d %d %d %d %d \n", nZ1, nZ0, nZ1+nZ0, n0, n1, n0+n1);
      //j1=48610;
      //      Rprintf("we have data[48610:48630]= %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", data[j1], data[j1+1], data[j1+2], data[j1+3], data[j1+4], data[j1+5], data[j1+6], data[j1+7], data[j1+8], data[j1+9], data[j1+10], data[j1+11], data[j1+12], data[j1+13], data[j1+14], data[j1+15], data[j1+16], data[j1+17], data[j1+18], data[j1+19]);
      //      Rprintf("and relative probX= %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", probX[j1], probX[j1+1], probX[j1+2], probX[j1+3], probX[j1+4], probX[j1+5], probX[j1+6], probX[j1+7], probX[j1+8], probX[j1+9], probX[j1+10], probX[j1+11], probX[j1+12], probX[j1+13], probX[j1+14], probX[j1+15], probX[j1+16], probX[j1+17], probX[j1+18], probX[j1+19]);
      //      Rprintf("and relative stateX= %d %d %d %d %d  %d %d %d %d %d %d %d %d %d %d  %d %d %d %d %d\n", stateX[j1], stateX[j1+1], stateX[j1+2], stateX[j1+3], stateX[j1+4], stateX[j1+5], stateX[j1+6], stateX[j1+7], stateX[j1+8], stateX[j1+9], stateX[j1+10], stateX[j1+11], stateX[j1+12], stateX[j1+13], stateX[j1+14], stateX[j1+15], stateX[j1+16], stateX[j1+17], stateX[j1+18], stateX[j1+19]);

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
      //while(1)
	//{
	  for(i=0;i<I;i++)
	    { 
	      a=AQ[i]+n[i][1];
 	      b=BQ[i]+n[i][0];
	      q[i]=rbeta(a,b);
	      //Rprintf("sample q %lf %lf %lf\n", a, b, q[i]);
	    }
	  //pi=stateX[0]*q[0]/(q[0]+1-q[1])+(1-stateX[0])*(1-q[1])/(q[0]+1-q[1]);
	  //if(pi>u1&q[0]>gzero&q[1]>gzero)//rejection sampling step
	  //if(q[0]>gzero&q[1]>gzero)
	  //  {
	   //   break;
	   // }
	//}
      //printf("q= %lf %lf\n", q[0],q[1]);
      
      // 2.2 sample inner mixture portion pi
      a=Api+nZ1;
      b=Bpi+nZ0;
      pi=rbeta(a,b);
      //printf("pi= %lf\n", pi);
      //pi=1;//no-zero inflated model
      //printf("q and pi= %lf %lf %lf\n", q[0],q[1],pi);

      // 2.3. sample lambda0 and mu1,phi1 when method="ZIP&NB"
      if (*met==0)
	{
	  // 2.3.1. sample lambda0
	  a=Alambda[0]+sumZ1;
	  b=1/(Blambda[0]+nZ1);
	  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda[0], sumZ1, Blambda[0], nZ1);
	  lambda[0]=rgamma(a,b);

	  // 2.3.2. sample mu1 from truncated norm(mu1, var_mu1, 0, infty)
	  //while(1)
	    //{
	      //propmu1=rnorm(mu[1],sdmu1);
	      //if (propmu1>0)
		//{
		  //break;
		//}
	    //}
	  // we use X~ F(X) then X is from F^-1(U), where U is a rv from U(0,1). Note CDF of truncated norm is (PHI((x-mu)/sigma)-PHI(-mu/sigma))/PHI(mu/sigma)
	  U1=runif(0,1);
	  temp=U1*pnorm(mu[1]/sdmu1, 0, 1, 1, 0)+pnorm(-mu[1]/sdmu1, 0, 1, 1, 0);
	  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu1+mu[1];	  
	  Cmu1=pnorm(mu[1]/sdmu1,0,1,1,1)-pnorm(propmu1/sdmu1,0,1,1,1);
	  logAmu1=sum1*(log(propmu1/(propmu1+phi[1]))-log(mu[1]/(mu[1]+phi[1])))+phi[1]*n1*(log(mu[1]+phi[1])-log(propmu1+phi[1]))+(Amu[1]-1.0)*log(propmu1/mu[1])-Bmu[1]*(propmu1-mu[1])+Cmu1;
	  U=runif(0,1);
	  if(log(U)<logAmu1)
	    {
	      mu[1]=propmu1;
	      rate[0]=rate[0]+1.0;
	    }
	  // 2.3.3. sample phi1
	  //while(1)
	    //{
	     // propphi1=rnorm(phi[1],sdphi1);
	      //if (propphi1>0)
		//{
		 // break;
		//}
	    //}
	  U1=runif(0,1);
	  temp=U1*pnorm(phi[1]/sdphi1, 0, 1, 1, 0)+pnorm(-phi[1]/sdphi1, 0, 1, 1, 0);
	  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi1+phi[1];	
	  Cphi1=pnorm(phi[1]/sdphi1,0,1,1,1)-pnorm(propphi1/sdphi1,0,1,1,1);
	  sumgammaphi1=0;
	  for (s=0;s<S;s++)
	    {
	      if ((stateX[s]==1)&(data[s]>0))
		    {
		       for (j2=0;j2<data[s];j2++)
		      {
		        sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi[1]+j2);
		      }
		    }
	    }
	  logAphi1=sumgammaphi1+sum1*(log(mu[1]+phi[1])-log(mu[1]+propphi1))+n1*(propphi1*log(propphi1/(mu[1]+propphi1))-phi[1]*log(phi[1]/(mu[1]+phi[1])))+(Aphi[1]-1.0)*log(propphi1/phi[1])-Bphi[1]*(propphi1-phi[1])+Cphi1;
	  U=runif(0,1);
	  if(log(U)<logAphi1)
	    {
	      phi[1]=propphi1;
	      rate[1]=rate[1]+1.0;
	    }
	}

      // 2.4. sample lambda when method="poisson"(*met=1)
      if (*met==1)
	{
	  a=Alambda[0]+sumZ1;
	  b=1/(Blambda[0]+nZ1);
	  //Rprintf("sample lambda0 %lf %lf %lf %d\n", Alambda[0], sumZ1, Blambda[0], nZ1);
	  lambda[0]=rgamma(a,b);
	  a=Alambda[1]+sum1;
	  b=1/(Blambda[1]+n1);
	  lambda[1]=rgamma(a,b);
	  //Rprintf("sample lambda1 %lf %lf %lf %d\n", Alambda[1], sum1, Blambda[1], n1);
	  if(lambda[0]>lambda[1])
	    {
	      temp=lambda[1];
	      lambda[1]=lambda[0];
	      lambda[0]=temp;
	      //Rprintf("lambda= %lf %lf\n", lambda[0], lambda[1]);
	    }
	  //Rprintf("lambda= %lf %lf\n", lambda[0], lambda[1]);
	}
      
      // 2.5. sample mu and phi when method="NB" (*met=2)
      if (*met==2)   
	{
	  // 2.5.1.Sample proposal mu from truncated normal distributions (postive tail)
	  //while(1)
	    //{
	      //propmu0=rnorm(mu[0],sdmu0);
	      //if(propmu0>0)
		//{
		  //break;
		//}
	    //}
	  //while(1)
	    //{
	      //propmu1=rnorm(mu[1],sdmu1);
	      //if (propmu1>0)
		//{
		  //break;
		//}
	    //}
	  U1=runif(0,1);
	  temp=U1*pnorm(mu[0]/sdmu0, 0, 1, 1, 0)+pnorm(-mu[0]/sdmu0, 0, 1, 1, 0);
	  propmu0=qnorm(temp, 0, 1, 1, 0)*sdmu0+mu[0];	
	  U1=runif(0,1);
	  temp=U1*pnorm(mu[1]/sdmu1, 0, 1, 1, 0)+pnorm(-mu[1]/sdmu1, 0, 1, 1, 0);
	  propmu1=qnorm(temp, 0, 1, 1, 0)*sdmu1+mu[1];	

	  // 2.5.2. Calculate the acceptance ratio for mu
	  Cmu0=pnorm(mu[0]/sdmu0,0,1,1,1)-pnorm(propmu0/sdmu0,0,1,1,1);// log(PHI(mu^t/sd_mu))-log(PHI(mu'/sd_mu))##using R function pnorm(x, mu, phi, lowtailistrue=1, logistrue=1)
	  Cmu1=pnorm(mu[1]/sdmu1,0,1,1,1)-pnorm(propmu1/sdmu1,0,1,1,1);
	  logAmu0=sumZ1*(log(propmu0/(propmu0+phi[0]))-log(mu[0]/(mu[0]+phi[0])))+phi[0]*nZ1*(log(mu[0]+phi[0])-log(propmu0+phi[0]))+(Amu[0]-1.0)*log(propmu0/mu[0])-Bmu[0]*(propmu0-mu[0])+Cmu0;
	  logAmu1=sum1*(log(propmu1/(propmu1+phi[1]))-log(mu[1]/(mu[1]+phi[1])))+phi[1]*n1*(log(mu[1]+phi[1])-log(propmu1+phi[1]))+(Amu[1]-1.0)*log(propmu1/mu[1])-Bmu[1]*(propmu1-mu[1])+Cmu1;
	  // 2.5.3. sample mu
	  U=runif(0,1);
	  if(log(U)<logAmu0)
	    {
	      mu[0]=propmu0;
	      rate[2]=rate[2]+1.0;
	    }
	  U=runif(0,1);
	  if(log(U)<logAmu1)
	    {
	      mu[1]=propmu1;
	      rate[0]=rate[0]+1.0;
	    }
	  //Rprintf("phi= %lf %lf\n", phi[0], phi[1]);
	
	  // 2.5.4. Sample proposal phi from truncated normal distribution (postive tail)
	  // while(1)
	    //{
	      //propphi0=rnorm(phi[0],sdphi0);
	      //if(propphi0>0)
		//{
		  //break;
		//}
	    //}
	  //while(1)
	    //{
	      //propphi1=rnorm(phi[1],sdphi1);
	      //if (propphi1>0)
		//{
		  //break;
		//}
	    //}
	  U1=runif(0,1);
	  temp=U1*pnorm(phi[0]/sdphi0, 0, 1, 1, 0)+pnorm(-phi[0]/sdphi0, 0, 1, 1, 0);
	  propphi0=qnorm(temp, 0, 1, 1, 0)*sdphi0+phi[0];	
	  U1=runif(0,1);
	  temp=U1*pnorm(phi[1]/sdphi1, 0, 1, 1, 0)+pnorm(-phi[1]/sdphi1, 0, 1, 1, 0);
	  propphi1=qnorm(temp, 0, 1, 1, 0)*sdphi1+phi[1];	

	  // 2.5.5. Calculate the acceptance ratio for phi
	  Cphi0=pnorm(phi[0]/sdphi0,0,1,1,1)-pnorm(propphi0/sdphi0,0,1,1,1);// log(PHI(phi^t/sd_phi))-log(PHI(phi'/sd_phi))##using R function pnorm(x, mu, phi, iflowtail=1, iflog=1)
	  Cphi1=pnorm(phi[1]/sdphi1,0,1,1,1)-pnorm(propphi1/sdphi1,0,1,1,1);
	  sumgammaphi0=0;
	  sumgammaphi1=0;
	  //double tempsum0=0;
	  //double tempsum1=0;
	  for (s=0;s<S;s++)
	    {
	      //if (stateX[s]==0&stateZ[s]==1)
	      //{
	  	  //sumgammaphi0=sumgammaphi0+lgammafn(data[s]+propphi0)+lgammafn(phi[0])-lgammafn(data[s]+phi[0])-lgammafn(propphi0);
	      //}
	      //if (stateX[s]==1)
	      //{
	  	  //sumgammaphi1=sumgammaphi1+lgammafn(data[s]+propphi1)+lgammafn(phi[1])-lgammafn(data[s]+phi[1])-lgammafn(propphi1);
	      //}
	      if ((stateX[s]==0)&(stateZ[s]==1)&(data[s]>0))
		    {
		      for (j2=0;j2<data[s];j2++)
		      {
		        sumgammaphi0=sumgammaphi0+log(propphi0+j2)-log(phi[0]+j2);
		      //tempsum0=tempsum0+log(propphi0+j2)-log(phi[0]+j2);
		      }
		    }
	      if ((stateX[s]==1)&(data[s]>0))
		    {
		      for (j2=0;j2<data[s];j2++)
		      {
		        sumgammaphi1=sumgammaphi1+log(propphi1+j2)-log(phi[1]+j2);
		      //tempsum1=tempsum1+log(propphi1+j2)-log(phi[1]+j2);
		      }
		    }		
	    }
	  //Rprintf("sumgamma0 and tempsum0 are %lf %lf\n", sumgammaphi0, tempsum0);	  
	  //Rprintf("sumgamma1 and tempsum1 are %lf %lf\n", sumgammaphi1, tempsum1);
	  logAphi0=sumgammaphi0+sumZ1*(log(mu[0]+phi[0])-log(mu[0]+propphi0))+nZ1*(propphi0*log(propphi0/(mu[0]+propphi0))-phi[0]*log(phi[0]/(mu[0]+phi[0])))+(Aphi[0]-1.0)*log(propphi0/phi[0])-Bphi[0]*(propphi0-phi[0])+Cphi0;
	  logAphi1=sumgammaphi1+sum1*(log(mu[1]+phi[1])-log(mu[1]+propphi1))+n1*(propphi1*log(propphi1/(mu[1]+propphi1))-phi[1]*log(phi[1]/(mu[1]+phi[1])))+(Aphi[1]-1.0)*log(propphi1/phi[1])-Bphi[1]*(propphi1-phi[1])+Cphi1;
	   // 2.5.6. sample phi
	  U=runif(0,1);
	  if(log(U)<logAphi0)
	    {
	      phi[0]=propphi0;
	      rate[3]=rate[3]+1.0;
	    }
	  U=runif(0,1);
	  if(log(U)<logAphi1)
	    {
	      phi[1]=propphi1;
	      rate[1]=rate[1]+1.0;
	    }
	  //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu[0], mu[1], phi[0], phi[1]);
	 
	  
	  if(mu[0]>mu[1])
	    {
	      //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu[0], mu[1], phi[0], phi[1]);
	      temp=phi[1];
	      phi[1]=phi[0];
	      phi[0]=temp;
	      temp=mu[1];
	      mu[1]=mu[0];
	      mu[0]=temp;
	    }
	  //Rprintf("mu,phi= %lf %lf %lf %lf\n", mu[0], mu[1], phi[0], phi[1]);
	}

      // iteration step 3. calculate likelihood function 
      logl=0.0;  
      pX1[1]=1/(1+exp(log(1-q[1])-log(q[0])));//stationary distribution
      pX1[0]=1-pX1[1];
      for (s=0; s<S; s++)
	{
	  logl=logl+log(((1-pi)*indexY[s]+pi*pY[0][s]*indexY[s]+pY[0][s]*(1-indexY[s]))*pX1[0]+pY[1][s]*pX1[1]);
	}
      
      // iteration step 4. pass the results back to R
      if ((z>=bN)&(z==j1*10+bN))
	{
	  *(es_pi+j1)=pi;
	  *(es_q0+j1)=q[0];
	  *(es_q1+j1)=q[1];
	  *(es_lambda1+j1)=lambda[1];
	  *(es_lambda0+j1)=lambda[0];
	  *(es_mu1+j1)=mu[1];
	  *(es_phi1+j1)=phi[1];
	  *(es_mu0+j1)=mu[0];
	  *(es_phi0+j1)=phi[0];
	  *(loglikeli+j1)=logl;  
	  j1++;
	  //Rprintf("%d %lf %lf %lf %lf %lf %lf %lf \n", z,q[1],q[0],mu[1],phi[1],pi,mu[0],phi[0]);
	  //Rprintf("%d %lf %lf %lf %lf %lf \n", z,q[1],q[0],lambda[1],pi,lambda[0]);
	}
      if (z==z1*jN)
	{
	  Rprintf("%d % ",z1*10);
	  z1++;	    
	}
    }  
  Rprintf("\n");
  PutRNGstate();
  double iter;
  iter=*N;
  //Rprintf("%lf %lf %lf %lf \n",rate[0]/iter,rate[1]/iter, rate[2]/iter, rate[3]/iter);// acceptance rate
  *(acrate+0)=rate[0]/iter;
  *(acrate+1)=rate[1]/iter;
  *(acrate+2)=rate[2]/iter;
  *(acrate+3)=rate[3]/iter;
 
  for (i=0; i<S; i++)
    {
      *(PP+i)=sumprobX[i]/iter;
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


