#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <complex.h>
#include "./solvers.h"
extern double ranf(void);
int main(int argc, char **argv)
{
  complex double csol[4]; 
  complex double c[5];
  int numtrials, its, caso=1;
  srand48(4242);
  numtrials=100000000;
  if (argc > 1)
    caso=atoi(argv[1]);
  if (argc > 2)
    numtrials = atoi(argv[2]);
  printf("Timing Test of ");
  switch (caso)
    {
    case 1:
      printf("OQS\n");
      break;
    default:
      printf("nothing\n");
      break;
    }
  for (its=0; its < numtrials; its++)
    {
      c[4]=1.0;
      c[3]=ranf()-0.5+I*(ranf()-0.5);
      c[2]=ranf()-0.5+I*(ranf()-0.5);
      c[1]=ranf()-0.5+I*(ranf()-0.5);
      c[0]=ranf()-0.5+I*(ranf()-0.5);
      switch (caso)
	{
	case 1:
	  oqs_quartic_solver_cmplx(c, csol);     
	  break;
	default:
	 // do nothing
	 break;
	}
    }
  exit(0);
}
