#include<stdlib.h>
#include<math.h>
#include <time.h>
#include<stdio.h>
#define Sqr(x) ((x)*(x))
#include <complex.h>
#include "./solvers.h"
#ifndef CMPLX
#define CMPLX(x,y) ((x)+I*(y))
#endif
#ifndef CMPLXL
#define CMPLXL(x,y) ((x)+I*(y))
#endif
extern int perm[24][4]; 
extern long double maxfact;
extern void sort_sol_opt(complex long double *sol, complex long double *exsol);
extern void sort_sol_optl(complex long double *sol, complex long double *exsol);
extern void oqs_solve_quadratic(long double a, long double b, complex long double roots[2]);
extern void oqs_quartic_solver_d20(long double coeff[5], complex long double roots[4], int chk);
void oqs_solve_quadratic_l(long double a, long double b, long double roots[2])
{ 
  long double div,diskr,zmax,zmin;
  diskr=a*a-4*b;   
  if(a>=0.0)
    div=-a-sqrt(diskr);
  else
    div=-a+sqrt(diskr);

  zmax=div/2;

  if(zmax==0.0)
    zmin=0.0;
  else
    zmin=b/zmax;

  roots[0]=zmax;
  roots[1]=zmin;
}
#ifdef OQS_LONG_DBL
long long double print_accuracy_atl(char *str, complex long long double *csol, complex long double exsol[4])
{
  /* we follow FLocke here */
  int k1;
  long long double relerr, relerrmax;
  for (k1=0; k1 < 4; k1++)
    {
      relerr=cabsl((csol[k1] - exsol[k1])/exsol[k1]); 
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=cabsl((csol[k1] - exsol[k1])/exsol[k1]); 
        }
    }
  printf("[%s] relative accuracy=%.16LG\n", str, relerrmax);
  return relerrmax;
}
#endif

long double print_accuracy_at(char *str, long double complex *csol, complex long double exsol[4])
{
  /* we follow FLocke here */
  int k1;
  long double relerr, relerrmax;
  for (k1=0; k1 < 4; k1++)
    {
      relerr=cabsl((csol[k1] - exsol[k1])/exsol[k1]); 
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=cabsl((csol[k1] - exsol[k1])/exsol[k1]); 
        }
    }
  printf("[%s] relative accuracy=%.16G\n", str, relerrmax);
  return relerrmax;
}

void print_roots(char *str, complex long double x1c, complex long double x2c,
                 complex long double x3c, complex long double x4c)
{
  printf("%s\n", str);
  printf("root #1= %.36LG+I*(%.36LG)\n", creall(x1c), cimagl(x1c));
  printf("root #2= %.36LG+I*(%.36LG)\n", creall(x2c), cimagl(x2c));
  printf("root #3= %.36LG+I*(%.36LG)\n", creall(x3c), cimagl(x3c));
  printf("root #4= %.36LG+I*(%.36LG)\n", creall(x4c), cimagl(x4c));
}
int main(int argc, char **argv)
{
  complex long double x1c, x2c, x3c, x4c; 
  complex long double csol[4];
  complex long double csolREF[4];
  long double c[5], A, B, err;
  long int i;
  int k1;
#ifdef OQS_LONG_DBL
  long long double cl[5];
  complex long long double csoll[4];
#endif

  srand48(time(NULL));
  x1c=x2c=x3c=x4c=0.0;
  int fine = 0;
  complex long double qr[2];
  long double tmpfact=1.0, epsfact = 1.0;
  int maxtrials = (argc > 1)?atoi(argv[1]):100000;
  for (i=0; i < maxtrials && !fine; i++)
    {
      A = ((drand48()-0.5)*1);
      B = ((drand48()-0.5)*1);
      //printf("CASE 26\n");

      c[4]=1.0;
      c[3]=0.0;
      c[2]=A;
      c[1]=0.0;
      c[0]=B;

#if 0
      oqs_solve_quadratic(A, B, qr);
     //printf("qr = %15G %.15G\n", qr[0], qr[1]);
      //printf("A=%.16G B=%.16G A*A - 4.0*B=%.16G\n", A, B, A*A-4*B);
      csolREF[0]=csqrt(qr[0]);
      csolREF[1]=-csqrt(qr[0]);
      csolREF[2]=csqrt(qr[1]);
      csolREF[3]=-csqrt(qr[1]);
#else
      oqs_quartic_solver_d20(c, csolREF, 1);     
      if (maxfact!=1.0)
        tmpfact = maxfact;
 #endif
      oqs_quartic_solver_d20(c, csol, 0);     
   
      sort_sol_optl(csol, csolREF);

      err = 0.0;
      for (k1=0; k1 < 4; k1++)
        {
          err += cabsl(csolREF[k1])!=0?(cabsl(csol[k1]-csolREF[k1])/cabsl(csolREF[k1])):cabsl(csol[k1]-csolREF[k1]);
        }
      if (err > 1)
        {
          if (tmpfact > epsfact)
            epsfact = tmpfact;
#if 0
          printf("BAD POLYNOMIAL\n");
          printf("x^4 + (%.16LG)*x^3 + (%.16LG)*x^2 + (%.1L6G)*x + %.16LG\n", c[3], c[2], c[1], c[0]);
          print_roots("found roots", csol[0], csol[1], csol[2], csol[3]);
          print_roots("exact roots", csolREF[0], csolREF[1], csolREF[2], csolREF[3]);
          printf("err=%16LG\n", err);
          printf("maxfact=%.16LG\n", epsfact);
#endif
          fine = 0;
          //exit(-1);
        }
    }
  printf("maxfact=%.16LG\n", epsfact);
  exit(-1);
}
