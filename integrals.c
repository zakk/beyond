#include "integrals.h"
#include "mathematica.h"
#include "cubature.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern struct mathlink_env_t *met;

extern unsigned maxEval;
extern double reqAbsError;
extern double reqRelError;

double mathlink_Ij(double mu,double delta,double T)
{
	const char *res;
	double ret;
	
	res=mathlink_eval(met,"NIntegrate[fIj[k,%f,%f,%f],{k,0,Infinity}]",mu,delta,T);
	ret=atof(res);
	mathlink_disown_string(met,res);

	return ret;
}

double mathlink_Ik(double mu,double delta,double T)
{
	const char *res;
	double ret;
	
	res=mathlink_eval(met,"NIntegrate[fIj[k,%f,%f,%f],{k,0,Infinity}]",mu,delta,T);
	ret=atof(res);
	mathlink_disown_string(met,res);

	return ret;
}

double mathlink_Ijprime(double mu,double delta,double T)
{
	const char *res;
	double ret;
	
	res=mathlink_eval(met,"NIntegrate[D[fIj[k,x,%f,%f],x] /. {x -> %f},{k,0,Infinity}]",delta,T,mu);
	ret=atof(res);
	mathlink_disown_string(met,res);

	return ret;
}

double mathlink_Ikprime(double mu,double delta,double T)
{
	const char *res;
	double ret;
	
	res=mathlink_eval(met,"NIntegrate[D[fIk[k,x,%f,%f],x] /. {x -> %f},{k,0,Infinity}]",delta,T,mu);
	ret=atof(res);
	mathlink_disown_string(met,res);

	return ret;
}

struct generic_ctx_t
{
	double mu,delta,T;
	double (*func)(double,double,double,double);
};

void fIgeneric(unsigned ndim,const double *in, void *fdata,
         unsigned fdim,double *fval)
{
	double k,y=in[0];
	double mu,delta,T;
	double ret;
	
	struct generic_ctx_t *ctx=(struct generic_ctx_t *)(fdata);

	mu=ctx->mu;
	delta=ctx->delta;
	T=ctx->T;

	k=(1.0f/y)-1;
	ret=ctx->func(k,mu,delta,T);

	fval[0]=ret/(y*y);
}

double Igeneric(double (*func)(double,double,double,double),double mu,double delta,double T)
{
	void *fdata;
	double xmin=0.0f;
	double xmax=1.0f;
	double ret,err;

	struct generic_ctx_t ctx;

	ctx.mu=mu;
	ctx.delta=delta;
	ctx.T=T;
	ctx.func=func;

	fdata=&ctx;
	
	adapt_integrate(1,fIgeneric,fdata,1,&xmin,&xmax,maxEval,reqAbsError,reqRelError,&ret,&err);

	return ret;
}

#define Power(x,y)	pow((double)(x),(double)(y))
#define Sqrt(z)		sqrt((double)(z))
#define Sinh(z)		sinh((double)(z))
#define Cosh(z)		cosh((double)(z))
#define Tanh(z)		tanh((double)(z))
#define Sech(z)		pow(cosh((double)(z)),-1.0f)
#define Csch(z)		pow(sinh((double)(z)),-1.0f)

double glue__fIj(double k,double mu,double delta,double T)
{
	double ret;

	ret=Power(k,2)*(-Power(k,2)/(2.*T*(1 + Cosh(Sqrt(Power(mu,2) + Power(Power(k,2) - mu,2))/(2.*T)))) + 
     	    (1 - ((Power(k,2) - mu)*Tanh(Sqrt(Power(mu,2) + Power(Power(k,2) - mu,2))/(2.*T)))/
	    Sqrt(Power(mu,2) + Power(Power(k,2) - mu,2)))/2.);

	return ret;
}

double glue__dfIj(double k,double mu,double delta,double T)
{
	double ret;

	ret=(Power(k,2)*(-((Power(k,2)*(Power(delta,2) + Power(Power(k,2) - mu,2))*
	            (Power(k,2) - mu)*Sinh(Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))/
	              (2.*T)))/
	          Power(1 + Cosh(Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))/(2.*T)),2))\
	        + T*Power(Sech(Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))/(2.*T)),2)*
	        (Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))*Power(Power(k,2) - mu,2) + 
	          T*Power(delta,2)*Sinh(Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))/T))))/
		(4.*Power(T,2)*Power(Power(delta,2) + Power(Power(k,2) - mu,2),1.5));

	return ret;
}

double glue__fIk(double k,double mu,double delta,double T)
{
	double ret;

	ret=Power(k,2)*(-(Power(Power(k,2) - mu,2)/
            (T*(Power(mu,2) + Power(Power(k,2) - mu,2))*
            (1 + Cosh(Sqrt(Power(mu,2) + Power(Power(k,2) - mu,2))/(2.*T))))) + 
            (Power(mu,2)*Tanh(Sqrt(Power(mu,2) + Power(Power(k,2) - mu,2))/(2.*T)))/
            Power(Power(mu,2) + Power(Power(k,2) - mu,2),1.5));

	return ret;
}

double glue__dfIk(double k,double mu,double delta,double T)
{
	double r1,r2,r3,r4,r5,r4a,r4b,r4c;

	r1=(2*(Power(k,2) - mu))/
	   (T*(Power(delta,2) + Power(Power(k,2) - mu,2))*
	(1 + Cosh(Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))/(2.*T))));

	r2=(2*Power(Power(k,2) - mu,3))/
	   (T*Power(Power(delta,2) + Power(Power(k,2) - mu,2),2)*
	     (1 + Cosh(Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))/(2.*T))));

	r3=(Power(delta,2)*(Power(k,2) - mu)*Power(Sech(Sqrt(Power(delta,2) + 
	          Power(Power(k,2) - mu,2))/(2.*T)),2))/
		(2.*T*Power(Power(delta,2) + Power(Power(k,2) - mu,2),2));

	r4a=Power(Power(k,2) - mu,3)*Sinh(Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))/(2.*T));
	r4b=2*Power(T,2)*Power(Power(delta,2) + Power(Power(k,2) - mu,2),1.5);
	r4c=Power(1 + Cosh(Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))/(2.*T)),2);

	r5=(3*Power(delta,2)*(Power(k,2) - mu)*Tanh(Sqrt(Power(delta,2) + Power(Power(k,2) - mu,2))/
	       (2.*T)))/Power(Power(delta,2) + Power(Power(k,2) - mu,2),2.5);

	r4=r4a/r4b/r4c;

	return k*k*(r1-r2-r3-r4+r5);
}

double cubature_Ij(double mu,double delta,double T)
{
	return Igeneric(glue__fIj,mu,delta,T);
}

double cubature_Ik(double mu,double delta,double T)
{
	return Igeneric(glue__fIk,mu,delta,T);
}

double cubature_Ijprime(double mu,double delta,double T)
{
	return Igeneric(glue__dfIj,mu,delta,T);
}

double cubature_Ikprime(double mu,double delta,double T)
{
	return Igeneric(glue__dfIk,mu,delta,T);
}
