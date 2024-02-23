#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<signal.h>
#include <complex.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#define Sqr(x) ((x)*(x))
#ifndef CMPLX
#define CMPLX(x,y) (x)+(y)*I
#endif
const long double cubic_rescal_fact = 3.488062113727083E+102; //= pow(DBL_MAX,1.0/3.0)/1.618034;
const long double quart_rescal_fact = 7.156344627944542E+76; // = pow(DBL_MAX,1.0/4.0)/1.618034;
const long double macheps =  5.421010862427522E-20; //2.2204460492503131E-16; // DBL_EPSILON
int oqs_check_always_d20 = 0;
long double oqs_max2(long double a, long double b)
{
  if (a >= b)
    return a;
  else
    return b;
}
long double oqs_max3(long double a, long double b, long double c)
{
  long double t;
  t = oqs_max2(a,b);
  return oqs_max2(t,c);
}

long double oqs_min2(long double a, long double b)
{
  if (a <= b)
    return a;
  else
    return b;
}

long double oqs_min3(long double a, long double b, long double c)
{
  long double t;
  t = oqs_min2(a,b);
  return oqs_min2(t,c);
}

void oqs_solve_cubic_analytic_depressed_handle_inf(long double b, long double c, long double *sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * where coefficients b and c are large (see sec. 2.2 in the manuscript) */ 
  long double Q, R, theta, A, B, QR, QRSQ, KK, sqrtQ, RQ;;
  const long double PI2=M_PI/2.0, TWOPI=2.0*M_PI;
  Q = -b/3.0;
  R = 0.5*c;
  if (R==0)
    {
      if (b <= 0)
        {
          *sol=sqrt(-b);
        }
      else
        {
          *sol=0;
        }
      return;
    }

  if (fabsl(Q) < fabsl(R))
    {
      QR=Q/R;
      QRSQ=QR*QR; 
      KK=1.0 - Q*QRSQ;
    }
  else
    {
      RQ = R/Q;
      KK = copysign(1.0,Q)*(RQ*RQ/Q-1.0);
    }

  if (KK < 0.0)
    {
      sqrtQ=sqrt(Q);
      theta = acos((R/fabsl(Q))/sqrtQ);
      if (theta < PI2) 
        *sol = -2.0*sqrtQ*cos(theta/3.0);
      else 
        *sol = -2.0*sqrtQ*cos((theta+TWOPI)/3.0);
    }
  else
    {
      if (fabsl(Q) < fabsl(R))
        A = -copysign(1.0,R)*cbrt(fabsl(R)*(1.0+sqrt(KK)));
      else
        {
          A = -copysign(1.0,R)*cbrt(fabsl(R)+sqrt(fabsl(Q))*fabsl(Q)*sqrt(KK));
        }
      if (A==0.0)
        B=0.0;
      else
        B = Q/A;
      *sol = A+B;
    }
}
void oqs_solve_cubic_analytic_depressed(long double b, long double c, long double *sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * (see sec. 2.2 in the manuscript) */ 
  long double Q, R, theta, Q3, R2, A, B, sqrtQ;
  Q = -b/3.0;
  R = 0.5*c;
  if (fabsl(Q) > 1E102 || fabsl(R) > 1E154)
    {
      oqs_solve_cubic_analytic_depressed_handle_inf(b, c, sol);
      return;
    }
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sqrtQ=-2.0*sqrt(Q);
      if (theta < M_PI/2) 
        *sol = sqrtQ*cos(theta/3.0);
      else 
        *sol = sqrtQ*cos((theta+2.0*M_PI)/3.0);
    }
  else
    {
      A = -copysign(1.0,R)*pow(fabsl(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
        B=0.0;
      else
        B = Q/A;
      *sol = A+B; /* this is always largest root even if A=B */
    }
}
long double oqs_calc_phi0(long double a, long double b, long double c, long double d, long double *phi0, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic 
   * in eq. (79) (see also the discussion in sec. 2.2 of the manuscript) */
  long double rmax, g,h,gg,hh,aq,bq,cq,dq,s,diskr;
  long double maxtt, xxx, gx, x, xold, f, fold, df, xsq;
  long double ggss, hhss, dqss, aqs, bqs, cqs, rfact, rfactsq; 
  int iter;
  diskr=9*a*a-24*b;                    
  /* eq. (87) */
  if(diskr > 0.0)
    { 
      diskr=sqrt(diskr);
      if(a > 0.0)
        s=-2*b/(3*a+diskr);                     
      else
        s=-2*b/(3*a-diskr);                      
    }
  else
    {      
      s=-a/4;                                    
    }
  /* eqs. (83) */
  aq=a+4*s;                                      
  bq=b+3*s*(a+2*s);                              
  cq=c+s*(2*b+s*(3*a+4*s));                      
  dq=d+s*(c+s*(b+s*(a+s)));                      
  gg=bq*bq/9;
  hh=aq*cq;     

  g=hh-4*dq-3*gg;                       /* eq. (85) */  
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq; /* eq. (86) */          
  oqs_solve_cubic_analytic_depressed(g, h, &rmax);
  if (isnan(rmax) || isinf(rmax))
    {
      oqs_solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
      if ((isnan(rmax) || isinf(rmax)) && scaled)
        {
          // try harder: rescale also the depressed cubic if quartic has been already rescaled
          rfact = cubic_rescal_fact; 
          rfactsq = rfact*rfact;
          ggss = gg/rfactsq;
          hhss = hh/rfactsq;
          dqss = dq/rfactsq;
          aqs = aq/rfact;
          bqs = bq/rfact;
          cqs = cq/rfact;
          ggss=bqs*bqs/9.0;
          hhss=aqs*cqs;   
          g=hhss-4.0*dqss-3.0*ggss;                       
          h=(8.0*dqss+hhss-2.0*ggss)*bqs/3-cqs*(cqs/rfact)-(dq/rfact)*aqs*aqs; 
          oqs_solve_cubic_analytic_depressed(g, h, &rmax);
          if (isnan(rmax) || isinf(rmax))
            {
              oqs_solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
            }
          rmax *= rfact;
        }
    }
  /* Newton-Raphson used to refine phi0 (see end of sec. 2.2 in the manuscript) */
  x = rmax;
  xsq=x*x;
  xxx=x*xsq;
  gx=g*x;
  f = x*(xsq + g) + h;
  if (fabsl(xxx) > fabsl(gx))
    maxtt = fabsl(xxx);
  else
    maxtt = fabsl(gx);
  if (fabsl(h) > maxtt)
    maxtt = fabsl(h);

  if (fabsl(f) > macheps*maxtt)
    {
      for (iter=0; iter < 8; iter++)
        {   
          df =  3.0*xsq + g;
          if (df==0)
            {
              break;
            }
          xold = x;
          x += -f/df;
          fold = f;
          xsq = x*x;
          f = x*(xsq + g) + h;
          if (f==0)
            {
              break;
            } 

          if (fabsl(f) >= fabsl(fold))
            {
              x = xold;
              break;
            }
        }
    }
  *phi0 = x;
  return f;
}
long double oqs_calc_err_ldlt(long double b, long double c, long double d, long double d2, long double l1, long double l2, long double l3)
{
  /* Eqs. (29) and (30) in the manuscript */
  long double sum;
  sum =  (b==0)?fabsl(d2 + l1*l1 + 2.0*l3):fabsl(((d2 + l1*l1 + 2.0*l3)-b)/b);
  sum += (c==0)?fabsl(2.0*d2*l2 + 2.0*l1*l3):fabsl(((2.0*d2*l2 + 2.0*l1*l3)-c)/c);
  sum += (d==0)?fabsl(d2*l2*l2 + l3*l3):fabsl(((d2*l2*l2 + l3*l3)-d)/d);
  return sum;
}
long double oqs_calc_err_abcd_cmplx(long double a, long double b, long double c, long double d, 
                               complex long double aq, complex long double bq, complex long double cq, complex long double dq)
{
  /* Eqs. (68) and (69) in the manuscript for complex alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq) */
  long double sum;
  sum = (d==0)?cabsl(bq*dq):cabsl((bq*dq-d)/d);
  sum += (c==0)?cabsl(bq*cq + aq*dq):cabsl(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?cabsl(bq + aq*cq + dq):cabsl(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?cabsl(aq + cq):cabsl(((aq + cq) - a)/a);
  return sum;
}
long double oqs_calc_err_d(long double errmin, long double d, long double bq, long double dq)
{
  /* Eqs. (68) and (69) in the manuscript for real alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq)*/
  return (d==0)?fabsl(bq*dq):fabsl((bq*dq-d)/d)+errmin;
}
long double oqs_calc_err_abcd(long double a, long double b, long double c, long double d, long double aq, long double bq, long double cq, long double dq)
{
  /* Eqs. (68) and (69) in the manuscript for real alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq)*/
  long double sum;
  sum = (d==0)?fabsl(bq*dq):fabsl((bq*dq-d)/d);
  sum += (c==0)?fabsl(bq*cq + aq*dq):fabsl(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?fabsl(bq + aq*cq + dq):fabsl(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?fabsl(aq + cq):fabsl(((aq + cq) - a)/a);
  return sum;
}
long double oqs_calc_err_abc(long double a, long double b, long double c, long double aq, long double bq, long double cq, long double dq)
{
  /* Eqs. (48)-(51) in the manuscript */
  long double sum;
  sum = (c==0)?fabsl(bq*cq + aq*dq):fabsl(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?fabsl(bq + aq*cq + dq):fabsl(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?fabsl(aq + cq):fabsl(((aq + cq) - a)/a);
  return sum;
}
#define NITERMAX 20
void oqs_NRabcd(long double a, long double b, long double c, long double d, long double *AQ, long double *BQ, long double *CQ, long double *DQ)
{
  /* Newton-Raphson described in sec. 2.3 of the manuscript for complex
   * coefficients a,b,c,d */
  int iter, k1, k2;
  long double x02, errfmin, errf, x[4], dx[4], det, Jinv[4][4], fvec[4];
  long double vr[4], errfold, errfa, fveca;
  long double xmin[4];
  x[0] = *AQ;
  x[1] = *BQ;
  x[2] = *CQ;
  x[3] = *DQ;
  vr[0] = d;
  vr[1] = c;
  vr[2] = b;
  vr[3] = a;
  fvec[0] = x[1]*x[3] - d;
  fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
  fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
  fvec[3] = x[0] + x[2] - a; 
  errf=0;
  errfa=0;
  for (k1=0; k1 < 4; k1++)
    {
      fveca = fabsl(fvec[k1]);
      errf += (vr[k1]==0)?fveca:fabsl(fveca/vr[k1]);
      errfa += fveca;
      xmin[k1] = x[k1];
    }
  errfmin = errfa;
  if (errfa==0)
    return;
  // on average 2 iterations are sufficient...
  for (iter = 0; iter < NITERMAX; iter++)
    {
      x02 = x[0]-x[2];
      det = x[1]*x[1] + x[1]*(-x[2]*x02 - 2.0*x[3]) + x[3]*(x[0]*x02 + x[3]);
      if (det==0.0)
        break;
      Jinv[0][0] = x02;
      Jinv[0][1] = x[3] - x[1];
      Jinv[0][2] = x[1]*x[2] - x[0]*x[3];
      Jinv[0][3] = -x[1]*Jinv[0][1] - x[0]*Jinv[0][2]; 
      Jinv[1][0] = x[0]*Jinv[0][0] + Jinv[0][1];
      Jinv[1][1] = -x[1]*Jinv[0][0];
      Jinv[1][2] = -x[1]*Jinv[0][1];   
      Jinv[1][3] = -x[1]*Jinv[0][2];
      Jinv[2][0] = -Jinv[0][0];
      Jinv[2][1] = -Jinv[0][1];
      Jinv[2][2] = -Jinv[0][2];
      Jinv[2][3] = Jinv[0][2]*x[2] + Jinv[0][1]*x[3];
      Jinv[3][0] = -x[2]*Jinv[0][0] - Jinv[0][1];
      Jinv[3][1] = Jinv[0][0]*x[3];
      Jinv[3][2] = x[3]*Jinv[0][1];
      Jinv[3][3] = x[3]*Jinv[0][2];
      for (k1=0; k1 < 4; k1++)
        {
          dx[k1] = 0;
          for (k2=0; k2 < 4; k2++)
            dx[k1] += Jinv[k1][k2]*fvec[k2];
        }
      for (k1=0; k1 < 4; k1++)
        {
          x[k1] += -dx[k1]/det;
        }
      fvec[0] = x[1]*x[3] - d;
      fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
      fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
      fvec[3] = x[0] + x[2] - a; 
      errfold = errf;
      errf=0;
      errfa = 0.0;
      for (k1=0; k1 < 4; k1++)
        {
          fveca=fabsl(fvec[k1]);
          errf += (vr[k1]==0)?fveca:fabsl(fveca/vr[k1]); 
          errfa += fveca;
        }

      if (errfa < errfmin)
        {
          errfmin = errfa;
          for (k1=0; k1 < 4; k1++)
            xmin[k1]=x[k1];
        }
      if (errfa==0)
        break;
      // stop if total relative has increased at least once but take solution with minimum absolute error
      if (errf >= errfold)
        break;
    }
  *AQ=xmin[0];
  *BQ=xmin[1];
  *CQ=xmin[2];
  *DQ=xmin[3];
}

void oqs_solve_quadratic(long double a, long double b, complex long double roots[2])
{ 
  long double div,sqrtd,diskr,zmax,zmin;
  diskr=a*a-4*b;   
  if(diskr>=0.0)
    {
      if(a>=0.0)
        div=-a-sqrt(diskr);
      else
        div=-a+sqrt(diskr);

      zmax=div/2;

      if(zmax==0.0)
        zmin=0.0;
      else
        zmin=b/zmax;

      roots[0]=CMPLX(zmax,0.0);
      roots[1]=CMPLX(zmin,0.0);
    } 
  else
    {   
      sqrtd = sqrt(-diskr);
      roots[0]=CMPLX(-a/2,sqrtd/2);
      roots[1]=CMPLX(-a/2,-sqrtd/2);      
    }   
}
long double maxfact=1.0;
void oqs_quartic_solver(long double coeff[5], complex long double roots[4])      
{
  /* USAGE:
   *
   * This routine calculates the roots of the quartic equation
   *
   * coeff[4]*x^4 + coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0] = 0
   * 
   * if coeff[4] != 0 
   *
   * the four roots will be stored in the complex array roots[] 
   *
   * */
  complex long double acx1, bcx1, ccx1, dcx1,acx=0.0+I*0.0,bcx=0.0+I*0.0,ccx,dcx,cdiskr,zx1,zx2,
          zxmax,zxmin, qroots[2];
  long double l2m[12], d2m[12], res[12], resmin, bl311, dml3l3, err0=0, err1=0, aq1, bq1, cq1, dq1; 
  long double a,b,c,d,phi0,aq,bq,cq,dq,d2,d3,l1,l2,l3, errmin, errv[3], aqv[3], cqv[3],gamma,del2,detM;
  int realcase[2], whichcase, k1, k, kmin, nsol;
  long double rfactsq, rfact=1.0, sqrtd3;

  if (coeff[4]==0.0)
    {
      printf("That's not a quartic!\n");
      return;
    }
  a=coeff[3]/coeff[4];
  b=coeff[2]/coeff[4];
  c=coeff[1]/coeff[4];
  d=coeff[0]/coeff[4];
  detM=oqs_calc_phi0(a,b,c,d,&phi0, 0);

  // simple polynomial rescaling
  if (isnan(phi0)||isinf(phi0))
    {
      rfact = quart_rescal_fact; 
      a /= rfact;
      rfactsq = rfact*rfact;
      b /= rfactsq;
      c /= rfactsq*rfact;
      d /= rfactsq*rfactsq;
      detM=oqs_calc_phi0(a,b,c,d,&phi0, 1);
    }
  l1=a/2;          /* eq. (16) */                                        
  l3=b/6+phi0/2;   /* eq. (18) */                                
  del2=c-a*l3;     /* defined just after eq. (27) */                             
  nsol=0;
  bl311 =2.*b/3.-phi0-l1*l1;   /* This is d2 as defined in eq. (20)*/ 
  dml3l3 = d-l3*l3;            /* dml3l3 is d3 as defined in eq. (15) with d2=0 */ 

  /* Three possible solutions for d2 and l2 (see eqs. (28) and discussion which follows) */
  long double d2min;
  if (bl311!=0.0)
    {
      d2m[nsol] = bl311;  
      l2m[nsol] = del2/(2.0*d2m[nsol]);   
      res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      d2min = d2m[nsol];
      nsol++;
    }
  if (del2!=0)
    {
      l2m[nsol]=2*dml3l3/del2;
      if (l2m[nsol]!=0)
        {
          d2m[nsol]=del2/(2*l2m[nsol]);
          res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
          if (fabsl(d2m[nsol]) < fabsl(d2min))
            d2min = d2m[nsol];
          nsol++;
        }

      d2m[nsol] = bl311;
      l2m[nsol] = 2.0*dml3l3/del2;
      res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      if (fabsl(bl311) < fabsl(d2min))
        d2min = bl311;
      nsol++;
    }

  if (nsol==0)
    {
      l2=d2=0.0;
    }
  else
    {
      /* we select the (d2,l2) pair which minimizes errors */
      for (k1=0; k1 < nsol; k1++)
        {
          if (k1==0 || res[k1] < resmin)
            {
              resmin = res[k1];
              kmin = k1;      
            }
        }
      d2 = d2m[kmin];
      l2 = l2m[kmin];
    }
  whichcase = 0; 
  if (d2 < 0.0) 
    {
      /* Case I eqs. (37)-(40) */
      gamma=sqrt(-d2);                               
      aq=l1+gamma;                                  
      bq=l3+gamma*l2;                              

      cq=l1-gamma;                                
      dq=l3-gamma*l2;                            
      if(fabsl(dq) < fabsl(bq))
        dq=d/bq;                                
      else if(fabsl(dq) > fabsl(bq))
        bq=d/dq;                               
      if (fabsl(aq) < fabsl(cq))
        {
          nsol=0;
          if (dq !=0)
            {
              aqv[nsol] = (c - bq*cq)/dq;    /* see eqs. (47) */
              errv[nsol]=oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
              nsol++;
            }
          if (cq != 0) 
            {
              aqv[nsol] = (b - dq - bq)/cq;  /* see eqs. (47) */
              errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
              nsol++;
            }
          aqv[nsol] = a - cq;                /* see eqs. (47) */
          errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
          nsol++;
          /* we select the value of aq (i.e. alpha1 in the manuscript) which minimizes errors */
          for (k=0; k < nsol; k++)
            {
              if (k==0 || errv[k] < errmin)
                {
                  kmin = k;
                  errmin = errv[k];
                }
            }
          aq = aqv[kmin];
        }
      else 
        {
          nsol = 0;
          if (bq != 0)
            { 
              cqv[nsol] = (c - aq*dq)/bq;              /* see eqs. (53) */
              errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
              nsol++;
            }
          if (aq != 0)
            {
              cqv[nsol] = (b - bq - dq)/aq;            /* see eqs. (53) */
              errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
              nsol++;
            }
          cqv[nsol] = a - aq;                          /* see eqs. (53) */
          errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
          nsol++;   
          /* we select the value of cq (i.e. alpha2 in the manuscript) which minimizes errors */
          for (k=0; k < nsol; k++)
            {
              if (k==0 || errv[k] < errmin)
                {
                  kmin = k;
                  errmin = errv[k];
                }
            }
          cq = cqv[kmin];
        }
      realcase[0]=1;
    }
  else if (d2 > 0)   
    {
      /* Case II eqs. (53)-(56) */
      gamma=sqrt(d2); 
      acx=CMPLX(l1,gamma);  
      bcx=CMPLX(l3,gamma*l2);
      ccx = conj(acx);
      dcx = conj(bcx);
      realcase[0] = 0; 
    }
  else 
    realcase[0] = -1; // d2=0
  /* Case III: d2 is 0 or approximately 0 (in this case check which solution is better) */
  // CRITICAL FIX 04/01/2024: 
  // if d2 == 0 then d3 != 0 and one has to consider 
  // the case 3 (d2 = 0) to find quartic roots (see pags. 9-10 of my ACM 2020).
  // To verify that d2==0 one can use Eq. (2) in my ACM remark, 
  // but in addition one can also check if d3 != 0 by using Eq. (15). 
  // To do this, we note that according to Eq.(15) 
  // d3 = det(M)/d2 = d - d2*l2^2 + l3^2, i.e.
  // det(M) = d2*d - d2^2*l2^2 + d2*l3^2
  // hence to consider d3 != 0 we require that
  // d3 > meps*min{abs(d2*d),abs(d2*d2*l2*l2),abs(l3*l3*d2)     
  
  // PREVIOUS CONDITION: if (realcase[0]==-1 || (fabsl(d2) <= macheps*oqs_max3(fabsl(2.*b/3.), fabsl(phi0), l1*l1))) 
  // FIX 29/12/2021: previous condition (see line above) was too stringent, hence I switched to criterion 2) in Ref. [28]
  //long double fdetM = fabsl(detM);
  //long double mepsd2 = fabsl(macheps*d2);
  //long double fd2 = fabsl(d2);
  //if (oqs_check_always_d20 || realcase[0]==-1 || (fabsl(d2) <= macheps*oqs_max3(fabsl(2.*b/3.), fabsl(phi0), l1*l1))) 
  //printf("detM=%.16G d - l3*l3=%.16G d2=%.16G detM/d2=%.16G\n", detM, d-l3*l3,d2,detM/d2);
  //long double d2alt = (2.0/3.0)*b - 0.25*a*a - phi0;  
  long double rhs= macheps*(fabsl(2.*b/3.)+fabsl(phi0)+l1*l1);
  //printf("rhs=%.16LG d2min=%.16LG\n", rhs, d2);
  //printf("sqrt(macheps)=%.16G\n", sqrt(macheps));
  maxfact=1.0;
  if (oqs_check_always_d20 || realcase[0]==-1 
      || fabsl(d2) <= macheps*(fabsl(2.*b/3.)+fabsl(phi0)+l1*l1))
    {
      maxfact = (d2!=0)?(fabsl(d2)/rhs):1.0;
      d3 = d - l3*l3;
      if (realcase[0]==1)
        err0 = oqs_calc_err_d(errmin, d, bq, dq);
      else if (realcase[0]==0)
        err0 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx, bcx, ccx, dcx);
      if (d3 <= 0)
        {
          realcase[1] = 1;
          sqrtd3 = sqrt(-d3);
          aq1 = l1;   
          bq1 = l3 + sqrtd3;
          cq1 = l1;
          dq1 = l3 - sqrtd3;
          if(fabsl(dq1) < fabsl(bq1))  
            dq1=d/bq1;                                        
          else if(fabsl(dq1) > fabsl(bq1))
            bq1=d/dq1;                                       
          err1 = oqs_calc_err_abcd(a, b, c, d, aq1, bq1, cq1, dq1); /* eq. (68) */
        }
      else /* complex */
        {
          realcase[1] = 0;
          acx1 = l1;
          bcx1 = l3 + I*sqrt(d3);
          ccx1 = l1;
          dcx1 = conj(bcx1);
          err1 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx1, bcx1, ccx1, dcx1); 
        }
      if (realcase[0]==-1 || err1 < err0)
        {
          whichcase=1; // d2 = 0
          if (realcase[1]==1)
            {
              aq = aq1;
              bq = bq1;
              cq = cq1;
              dq = dq1;
            }
          else
            {
              acx = acx1;
              bcx = bcx1;
              ccx = ccx1;
              dcx = dcx1;
            }
        }
    }
  if (realcase[whichcase]==1)
    {
      /* if alpha1, beta1, alpha2 and beta2 are real first refine 
       * the coefficient through a Newton-Raphson */
      oqs_NRabcd(a,b,c,d,&aq,&bq,&cq,&dq);      
      /* finally calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
      oqs_solve_quadratic(aq,bq,qroots);
      roots[0]=qroots[0];
      roots[1]=qroots[1];        
      oqs_solve_quadratic(cq,dq,qroots);
      roots[2]=qroots[0];
      roots[3]=qroots[1];
    }
  else
    {
      /* complex coefficients of p1 and p2 */
      if (whichcase==0) // d2!=0
        {
          cdiskr=acx*acx/4-bcx;               
          /* calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
          zx1=-acx/2+csqrt(cdiskr);
          zx2=-acx/2-csqrt(cdiskr);
          if(cabsl(zx1) > cabsl(zx2))
            zxmax=zx1;
          else
            zxmax=zx2;
          zxmin=bcx/zxmax;        
          roots[0]=zxmin;
          roots[1]=conj(zxmin);
          roots[2]=zxmax;
          roots[3]=conj(zxmax);
        }
      else // d2 ~ 0
        {
          /* never gets here! */
          cdiskr=csqrt(acx*acx-4.0*bcx);
          zx1 = -0.5*(acx+cdiskr);
          zx2 = -0.5*(acx-cdiskr);
          if (cabsl(zx1) > cabsl(zx2))
            zxmax = zx1;
          else
            zxmax = zx2;
          zxmin = bcx/zxmax;
          roots[0] = zxmax;
          roots[1] = zxmin;
          cdiskr=csqrt(ccx*ccx-4.0*dcx);
          zx1 = -0.5*(ccx+cdiskr);
          zx2 = -0.5*(ccx-cdiskr);
          if (cabsl(zx1) > cabsl(zx2))
            zxmax = zx1;
          else
            zxmax = zx2;
          zxmin = dcx/zxmax;
          roots[2]= zxmax;
          roots[3]= zxmin;
        }
    }
  if (rfact!=1.0)
    {
      for (k=0; k < 4; k++)
        roots[k] *= rfact;
    }
}
void oqs_quartic_solver_d20(long double coeff[5], complex long double roots[4], int chk)
{
  oqs_check_always_d20 = chk;
  oqs_quartic_solver(coeff, roots);
}
