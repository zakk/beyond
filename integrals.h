#ifndef INTEGRALS_H
#define	INTEGRALS_H

double mathlink_Ij(double mu,double delta,double T);
double mathlink_Ik(double mu,double delta,double T);
double mathlink_Ijprime(double mu,double delta,double T);
double mathlink_Ikprime(double mu,double delta,double T);

double cubature_Ij(double mu,double delta,double T);
double cubature_Ik(double mu,double delta,double T);
double cubature_Ijprime(double mu,double delta,double T);
double cubature_Ikprime(double mu,double delta,double T);

#ifdef MATHEMATICA_INTEGRALS
#define mathlink_Ij Ij
#define mathlink_Ik Ik
#define mathlink_Ijprime Ijprime
#define mathlink_Ikprime Ikprime
#else
#define cubature_Ij Ij
#define cubature_Ik Ik
#define cubature_Ijprime Ijprime
#define cubature_Ikprime Ikprime
#endif

#endif	/* INTEGRALS_H */
