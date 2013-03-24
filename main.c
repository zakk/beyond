#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#define MATHEMATICA_INTEGRALS

#include "mathematica.h"
#include "cubature.h"
#include "brent.h"
#include "integrals.h"

/*
	Questi parametri controllano la precisione dell'integrazione numerica degli integrali.
*/

unsigned maxEval=1024*128;
double reqAbsError=0;
double reqRelError=0;

double xi(double k,double mu,double delta)
{
	return k*k-mu;
}

double E(double k,double mu,double delta)
{
	return sqrt(pow(xi(k,mu,delta),2.0f)+delta*delta);
}

void fI1(unsigned ndim,const double *in, void *fdata,
         unsigned fdim,double *fval)
{
	double k,y=in[0];
	double mu=((double *)(fdata))[0];
	double delta=((double *)(fdata))[1];
	double T=((double *)(fdata))[2];
	double th,ret;

	k=(1.0f/y)-1;
	th=tanh(E(k,mu,delta)/(2.0f*T));
	ret=(k*k)*((th/E(k,mu,delta))-(1.0f/(k*k)));

	fval[0]=ret/(y*y);
}

double I1(double mu,double delta,double T)
{
	double fdata[3];
	double xmin=0.0f;
	double xmax=1.0f;
	double ret,err;

	fdata[0]=mu;
	fdata[1]=delta;
	fdata[2]=T;

	adapt_integrate(1,fI1,&fdata,1,&xmin,&xmax,maxEval,reqAbsError,reqRelError,&ret,&err);

	return ret;
}

void fI2(unsigned ndim,const double *in, void *fdata,
         unsigned fdim,double *fval)
{
	double k,y=in[0];
	double mu=((double *)(fdata))[0];
	double delta=((double *)(fdata))[1];
	double T=((double *)(fdata))[2];
	double th,ret;

	k=(1.0f/y)-1;
	th=tanh(E(k,mu,delta)/(2.0f*T));
	ret=(k*k)*(1.0f-(th*xi(k,mu,delta)/E(k,mu,delta)));

	fval[0]=ret/(y*y);
}

double I2beyond(double mu,double delta,double T);

double I2(double mu,double delta,double T)
{
	double fdata[3];
	double xmin=0.0f;
	double xmax=1.0f;
	double ret,err;

	fdata[0]=mu;
	fdata[1]=delta;
	fdata[2]=T;

	adapt_integrate(1,fI2,&fdata,1,&xmin,&xmax,maxEval,reqAbsError,reqRelError,&ret,&err);

	ret-=I2beyond(mu,delta,T);

	return ret;
}

double Ij(double mu,double delta,double T);
double Ik(double mu,double delta,double T);
double Ijprime(double mu,double delta,double T);
double Ikprime(double mu,double delta,double T);

double I2beyond(double mu,double delta,double T)
{
	double alpha,beta,aj,ak,ajprime,akprime;
	
	aj=Ij(mu,delta,T);
	ak=Ik(mu,delta,T);
	ajprime=Ijprime(mu,delta,T);
	akprime=Ikprime(mu,delta,T);

	alpha=pow(M_PI*T,4.0)/15.0f;
	beta=(ak*ajprime-aj*akprime)/(aj*ak);

	return alpha*beta;
}

void fIcond(unsigned ndim,const double *in, void *fdata,
         unsigned fdim,double *fval)
{
	double k,y=in[0];
	double mu=((double *)(fdata))[0];
	double delta=((double *)(fdata))[1];
	double T=((double *)(fdata))[2];
	double th,ret;

	k=(1.0f/y)-1;
	th=tanh(E(k,mu,delta)/(2.0f*T));
	ret=(k*k)*pow(th*delta/E(k,mu,delta),2.0)/4.0f;

	fval[0]=ret/(y*y);
}

double Icond(double mu,double delta,double T)
{
	double fdata[3];
	double xmin=0.0f;
	double xmax=1.0f;
	double ret,err;

	fdata[0]=mu;
	fdata[1]=delta;
	fdata[2]=T;

	adapt_integrate(1,fIcond,&fdata,1,&xmin,&xmax,maxEval,reqAbsError,reqRelError,&ret,&err);

	return ret;
}

double get_delta(double mu,double T)
{
	int status=0;
	double arg,value,machep;
	
	arg=0.0;
	machep=r8_epsilon();

	do
	{
		value=(3.0f/2.0f)*I2(mu,arg,T)-1.0f;
		zero_rc(1e-6,1e2,machep,&arg,&status,value);
	}
	while(status!=0);
	
	return arg;
}

double get_mu(double delta,double T)
{
	int status=0;
	double arg,value,machep;
	
	arg=0.0;
	machep=r8_epsilon();

	do
	{
		value=(3.0f/2.0f)*I2(arg,delta,T)-1.0f;
		zero_rc(-40.0f,1.0f,machep,&arg,&status,value);
	}
	while(status!=0);

	return arg;
}

double get_y(double mu,double delta,double T)
{	
	return -(2.0f/M_PI)*I1(mu,delta,T);
}

double get_phi(double mu,double delta,double T)
{
	return 3.0f*Icond(mu,delta,T);
}

struct mathlink_env_t *met;

/**
 * @brief Blah!;
 * @return nothing, it's void
 */

void cleanup_mathlink(void)
{
	fini_mathlink(met);
}

int main(void)
{
	double mu;
	double epsilon=10e-5;
	int cnt=0;

	FILE *out;
	char *outfile="jlevel1.dat";
	
#ifdef MATHEMATICA_INTEGRALS
	
	const char *fij="k^2*(-k^2/(2*t*(1 + Cosh[Sqrt[d^2 + (k^2 - m)^2]/(2*t)]))+"
		        "(1 - ((k^2 - m)*Tanh[Sqrt[d^2 + (k^2 - m)^2]/(2*t)])/Sqrt[d^2 + (k^2 - m)^2])/2)";

	const char *fik="k^2*(-((k^2 - m)^2/((d^2 + (k^2 - m)^2)*t*(1 + Cosh[Sqrt[d^2 + (k^2 - m)^2]/(2*t)]))) + "
		        "(d^2*Tanh[Sqrt[d^2 + (k^2 - m)^2]/(2*t)])/(d^2 + (k^2 - m)^2)^(3/2))";

	met=init_mathlink();
	atexit(cleanup_mathlink);

	mathlink_eval(met,"fIj[k_,m_,d_,t_] := %s",fij);
	mathlink_eval(met,"fIk[k_,m_,d_,t_] := %s",fik);

	fprintf(stderr,"Using Wolfram Mathematica for numerical integration.\n");
	
#endif
	
	if(!(out=fopen(outfile,"w+")))
	{
		printf("Couldn't open %s for writing.\n",outfile);
		return 0;
	}

	printf("Writing output to: %s\n",outfile);		

	for(mu=0.99;mu>=-5;mu-=0.5,cnt++)
	{	
		double cvals[40][2];
		double minT,maxT,deltaT;
		
		minT=epsilon;
		maxT=0.4+epsilon;
		deltaT=(maxT-minT)/((float)(40));

		for(int ticks=0;ticks<40;ticks++)
		{
			double T,delta,ij;
		
			T=minT+deltaT*ticks;
			delta=get_delta(mu,T);
			ij=Ij(mu,delta,T);

			cvals[ticks][0]=ij/(T*4*M_PI*0.45420f)-1;
			cvals[ticks][1]=T;

			printf("(%d/40 :: %d/12) ",ticks,cnt);
			fflush(stdout);
		}
	
		for(int ticks=1;ticks<40;ticks++)
		{
			if(cvals[ticks-1][0]*cvals[ticks][1]<0.0f)
			{
				double c1,c2;
				double T1,T2;
				double Tc,delta,y;
				
				c1=cvals[ticks-1][0];
				c2=cvals[ticks][0];
				T1=cvals[ticks-1][1];
				T2=cvals[ticks][1];

				Tc=T1+c1*(T2-T1)/(c2-c1);
				
				delta=get_delta(mu,Tc);
				y=get_y(mu,delta,Tc);
			
				fprintf(out,"%f %f %f %f\n",y,mu,delta,Tc);
				fflush(out);
			}
		}
	}

	if(out)
		fclose(out);
	
	return 0;
}

int main4(void)
{
	double ar[30][2]=
	{{0.14001, 1.}, {0.15001, 0.977989}, {0.16001, 0.99324}, {0.17001, 
  0.985887}, {0.18001, 1.}, {0.19001, 1.}, {0.20001, 
  0.999999}, {0.21001, 
  0.263496}, {0.22001, -0.245266}, {0.23001, -0.550679}, {0.24001, \
-0.783519}, {0.25001, -0.983729}, {0.26001, -1.16655}, {0.27001, \
-1.33891}, {0.28001, -1.50452}, {0.29001, -1.66557}};

	for(int i=0;i<30;i++)
	{
		double T=ar[i][0];
		double mu=ar[i][1];
		double delta=get_delta(mu,T);
		double y=get_y(mu,delta,T);
		
		printf("%f %f %f %f\n",y,mu,delta,T);
		fflush(stdout);
	}
	
	return 0;
}

int main3(void)
{
	for(double mu=-5.0;mu<=-0.95;mu+=0.1)
	{
		double T;
		double p;
		short first=1;
		
		for(T=0.0f;T<=0.5f;T+=0.02)
		{
			double delta=get_delta(mu,T);
			double y=get_y(mu,delta,T);
			double ij=Ij(mu,delta,T);
			double critical_value=ij/(T*4*M_PI*0.45420f);

			if((first==0)&&((critical_value-1)*(p-1)<0.0f))
			{
				printf("%f %f\n",y,T);
				fflush(stdout);
			}

			first=0;
			p=critical_value;
		}
	}
	
	return 0;
}


int main2(void)
{

#ifdef MATHEMATICA_INTEGRALS
	
	const char *fij="k^2*(-k^2/(2*t*(1 + Cosh[Sqrt[d^2 + (k^2 - m)^2]/(2*t)]))+"
		        "(1 - ((k^2 - m)*Tanh[Sqrt[d^2 + (k^2 - m)^2]/(2*t)])/Sqrt[d^2 + (k^2 - m)^2])/2)";

	const char *fik="k^2*(-((k^2 - m)^2/((d^2 + (k^2 - m)^2)*t*(1 + Cosh[Sqrt[d^2 + (k^2 - m)^2]/(2*t)]))) + "
		        "(d^2*Tanh[Sqrt[d^2 + (k^2 - m)^2]/(2*t)])/(d^2 + (k^2 - m)^2)^(3/2))";

	met=init_mathlink();
	atexit(cleanup_mathlink);

	mathlink_eval(met,"fIj[k_,m_,d_,t_] := %s",fij);
	mathlink_eval(met,"fIk[k_,m_,d_,t_] := %s",fik);

	fprintf(stderr,"Using Wolfram Mathematica for numerical integration.\n");
	
#endif

	FILE *out;
	char *outfile="jlevel0.dat";
	
	if(!(out=fopen(outfile,"w+")))
	{
		printf("Couldn't open %s for writing.\n",outfile);
		return 0;
	}

	printf("Writing output to: %s\n",outfile);
	
	int cnt=0;
	double mu,T;
	
	/*
	 *	ToDO: Tc binary search!
	 * 
	 */
	
	//for(int mucents=100;mucents>=0;mucents--)
	for(int mucents=0;mucents<=100;mucents++)
	{
		mu=0.95f-(0.95f+5.0f)/100.0f*mucents;
		
		double previous_critical_value,previous_T,previous_y;
		short first=1,exitloop=0;

		for(int Tcents=1;Tcents<=32;Tcents+=2)
		{
			T=Tcents*0.01f;

			double delta=get_delta(mu,T);
			double y=get_y(mu,delta,T);
			double ij=Ij(mu,delta,T);
			double critical_value=ij/(T*4*M_PI*0.45420f);

			if(first!=1)
			{
				double a,b;
				
				a=critical_value;
				b=previous_critical_value;

				if(((a-1.0f)*(b-1.0f))<0.0f)
				{
					double Tc=previous_T+(1.0f-b)*(T-previous_T)/(a-b);
					double yc=previous_y+(1.0f-b)*(y-previous_y)/(a-b);

					fprintf(out,"# T_C: %f %f\n",Tc,yc);
					fflush(out);
					exitloop=1;

				}
			}

			fprintf(out,"%f %f %f %f\n",y,mu,T,critical_value);
			fflush(out);

			previous_critical_value=critical_value;
			previous_T=T;
			previous_y=y;

			first=0;
			
			if(exitloop==1)
				break;
				
		}

		printf("%d/100 ",++cnt);
		fflush(stdout);
	}

	if(out)
		fclose(out);
	
	/*
	double delta,T;
	
	T=0.0f;
	for(delta=0.1f;delta<=2.0;delta+=0.1f)
	{
		double mu=get_mu(delta,T);
		double phi=get_phi(mu,delta,T);
		double y=get_y(mu,delta,T);

		printf("%f %f %f %f\n",y,delta,mu,phi);
		fflush(stdout);
	}
	*/
	
	/*
	double mu,T;
	
	T=0.0f;
	for(mu=0.95f;mu>=-5.0;mu-=0.1f)
	{
		double delta=get_delta(mu,T);
		double phi=get_phi(mu,delta,T);
		double y=get_y(mu,delta,T);

		printf("%f %f %f %f\n",y,delta,mu,phi);
		fflush(stdout);
	}
	*/
	
	/*
	double mu,T;

	for(mu=0.95f;mu>=-5.0;mu-=0.1f)
	{
		for(T=0.01;T<=1.0f;T+=0.01)
		{
			double delta=get_delta(mu,T);
			
			if(get_phi(mu,delta,T)<10e-3)
			{
				printf("%f %f\n",get_y(mu,delta,T),T);
				fflush(stdout);
				break;
			}
		}
	}
	*/

	/*
	double delta_threshold=10e-3;
	int tcents;

	for(tcents=0.0;tcents<=100;tcents++)
	{
		double T=tcents*0.01;
		double mu=get_mu(delta_threshold,T);
		double y=get_y(mu,delta_threshold,T);

		printf("%f %f\n",y,T);
		fflush(stdout);
	}
	*/

	/*
	double mu,T;
	
	T=0.2f;
	for(mu=0.95f;mu>=-5.0;mu-=0.1f)
	{
		double delta=get_delta(mu,T);
		double phi=get_phi(mu,delta,T);
		double y=get_y(mu,delta,T);

		printf("%f %f %f %f\n",y,delta,mu,phi);
		fflush(stdout);
	}
	*/
	/*	
	double delta,T;
	
	for(delta=0.1f;delta<=2.0;delta+=0.1f)
	{
		for(T=0;T<=1.0f;T+=0.1)
		{
			double mu=get_mu(delta,T);
			double y=get_y(mu,delta,T);
			
			printf("%f %f %f\n",y,T,Ij(mu,delta,T));
			fflush(stdout);
		}
	}
*/

	return 0;
}
