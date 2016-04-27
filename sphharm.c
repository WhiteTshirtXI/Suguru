#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "c99.h"
#include "name.h"
//for gsl library
//#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_sf_bessel.h>
#include <gsl_sf_exp.h>
#include <gsl_sf_gamma.h>
#include <gsl_sf_log.h>
#include <gsl_sf_pow_int.h>
#include <gsl_sf_legendre.h>
//for gsl library


#define fcal_legendre FORTRAN_NAME(cal_legendre,CAL_LEGENDRE)
#define fcal_legendre_array FORTRAN_NAME(cal_legendre_array,CAL_LEGENDRE_ARRAY)
#define fcal_legendre_deri1_array FORTRAN_NAME(cal_legendre_deri1_array,CAL_LEGENDRE_DERI1_ARRAY)
#define fcal_legendre_deri2_array FORTRAN_NAME(cal_legendre_deri2_array,CAL_LEGENDRE_DERI2_ARRAY)
#define fcal_legendre_deri4_array FORTRAN_NAME(cal_legendre_deri4_array,CAL_LEGENDRE_DERI4_ARRAY)
#define fcal_legendre_deri_array_asso FORTRAN_NAME(cal_legendre_deri_array_asso,CAL_LEGENDRE_DERI_ARRAY_ASSO)
//interface to use the gsl package to calculate Associated Legendre Polynomials and Spherical Harmonics
//NOTE that this Associated Legendre Polynomials and Spherical Harmonics has different definition with the one in Zhao's paper. Zhao's definition is from the library SPHEREPACK.
//the definition here is shown here http://www.gnu.org/software/gsl/manual/html_node/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html
#define fcal_deri1_array FORTRAN_NAME(cal_deri1_array,CAL_DERI1_ARRAY)


void fcal_legendre(double * legcoeff, int * l, int * m, double * x)
{
  *legcoeff = gsl_sf_legendre_sphPlm(*l,*m,*x)*sqrt(2*M_PI)/pow(-1,*m);
}
// interface to use the gsl package to calculate Associated Legendre Polynomials and Spherical Harmonics

void fcal_legendre_array(double * pbar, int * lmax, int * m, double * x)
{
  int ell;
  int cm = *m; int clmax = *lmax;
  double cx = *x;
  double coef=sqrt(2*M_PI)/pow(-1,cm);
  //  printf("%d %d %f\n",clmax,cm,cx);
  gsl_sf_legendre_sphPlm_array(clmax,cm,cx,pbar);
  //  printf("finish here\n");
  for(ell = cm ;ell <= clmax; ell++)
  pbar[ell - cm]      *= coef;
  //  printf("finish here22222\n");
}



void fcal_legendre_deri1_array(double * pbar, double * pbar_deri, int * lmax, int * m, double * x)
{
  int ell;
  int cm = *m; int clmax = *lmax;
  double cx = *x;
  double coef=sqrt(2*M_PI)/pow(-1,cm);
  //legcoeff has a size of *lmax
  // the size of pbar, pbar_deri and pbar_deri2 is *lmax - *m +1
  gsl_sf_legendre_sphPlm_deriv_array(clmax,cm,cx,pbar,pbar_deri);
  for(ell = cm ;ell <= clmax; ell++)
    {
      pbar[ell - cm]      *= coef;
      pbar_deri[ell - cm] *= coef;
    }
}

void fcal_legendre_deri2_array(double * pbar, double * pbar_deri, double * pbar_deri2, int * lmax, int * m, double * x)
{
  int ell;
  int cm = *m; int clmax = *lmax;
  double cx = *x;
  double coef=sqrt(2*M_PI)/pow(-1,cm);
  double coef2 = 1./(pow(cx,2)-1);
  //legcoeff has a size of *lmax
  // the size of pbar, pbar_deri and pbar_deri2 is *lmax - *m +1
  gsl_sf_legendre_sphPlm_deriv_array(clmax,cm,cx,pbar,pbar_deri);
  for(ell = cm ;ell <= clmax; ell++)
    {
      pbar[ell - cm]      *= coef;
      pbar_deri[ell - cm] *= coef;
    }

  pbar_deri2[0] = coef2*(cm*pbar[0]+(cm-2)*cx*pbar_deri[0]);
  
  for (ell = cm+1; ell <= clmax; ell++)
    {
      pbar_deri2[ell-cm]=coef2*(ell*pbar[ell-cm] + (ell-2)*cx*pbar_deri[ell-cm] - sqrt((2*ell+1)*(pow(ell,2)-pow(cm,2))/(2*ell-1))*pbar_deri[ell-cm-1]);
    }


  
  //////test of the recursion relation of the normalized spherical harmonic//
  /* ell = *m; */
  /* deri_test[ell - *m] = ell*(*x)/(pow(*x,2)-1) * pbar[ell - *m]; */
  /* for (ell = (*m)+1; ell <= *lmax; ell++) */
  /*   { */
  /*     deri_test[ell - *m] = (ell*(*x) * pbar[ell - *m] - sqrt((2*ell+1)*(pow(ell,2)-pow(*m,2))/(2*ell-1)) * pbar[ell - *m - 1])/(pow(*x,2)-1); */
  /*   } */

  /* for(ell = (*m) ;ell <= *lmax;ell++) */
  /*   { */
  /*     printf("%f \n", pbar_deri[ell - *m] - deri_test[ell - *m]); */
  /*   } */
  //////test of the recursion relation of the normalized spherical harmonic//
}

void fcal_legendre_deri4_array(double * pbar, double * pbar_deri1, double * pbar_deri2, double * pbar_deri3, double * pbar_deri4, int * lmax, int * m, double * x)
{
  int ell;
  int cm = *m; int clmax = *lmax;
  double cx = *x;
  double coef=sqrt(2*M_PI)/pow(-1,cm);
  double coef2 = 1./(pow(cx,2)-1);
  double coef3;
  int ind;
  //legcoeff has a size of *lmax
  // the size of pbar, pbar_deri1 and pbar_deri2 is *lmax - *m +1
  gsl_sf_legendre_sphPlm_deriv_array(clmax,cm,cx,pbar,pbar_deri1);
  for(ell = cm ;ell <= clmax; ell++)
    {
      pbar[ell - cm]      *= coef;
      pbar_deri1[ell - cm] *= coef;
    }

  pbar_deri2[0] = coef2*((cm    )*pbar[0]       + (cm-2)*cx*pbar_deri1[0]);
  pbar_deri3[0] = coef2*((2*cm-2)*pbar_deri1[0] + (cm-4)*cx*pbar_deri2[0]);
  pbar_deri4[0] = coef2*((3*cm-6)*pbar_deri2[0] + (cm-6)*cx*pbar_deri3[0]);

  for (ell = cm+1; ell <= clmax; ell++)
    {
      ind = ell - cm;
      coef3 = sqrt((2*ell+1)*(pow(ell,2)-pow(cm,2))/(2*ell-1)); 

      pbar_deri2[ind]=coef2*(ell*pbar[ind] + (ell-2)*cx*pbar_deri1[ind] - coef3*pbar_deri1[ind-1]);

      pbar_deri3[ind]=coef2*((2*ell-2)*pbar_deri1[ind] + (ell-4)*cx*pbar_deri2[ind] - coef3*pbar_deri2[ind-1]);

      pbar_deri4[ind]=coef2*((3*ell-6)*pbar_deri2[ind] + (ell-6)*cx*pbar_deri3[ind] - coef3*pbar_deri3[ind-1]);
    }


}

void fcal_legendre_deri_array_asso(int * lmax, int *m, double *x)
{
  int clmax = *lmax;
  int cm    = *m;
  int cx    = *x;
  double P[clmax-cm+1];
  double P_deri[clmax-cm+1], P_deri_my[clmax-cm+1];
  int ell;
  gsl_sf_legendre_Plm_deriv_array(clmax,cm,cx,P,P_deri);

  P_deri_my[0] = cm*cx/(pow(cx,2)-1)*P[0];
  for (ell = cm+1;ell <= clmax; ell ++)
    {
      P_deri_my[ell-cm] = 1./(pow(cx,2)-1)*(ell*cx*P[ell-cm] - (ell+cm)*P[ell-cm-1]);
    }
  for (ell = cm; ell <= clmax; ell++)
    {
      printf("%f %f\n",P_deri_my[ell], P_deri[ell]);
    }
}
 
void fcal_deri1(double * pbar, int l, int * m, double * x)
{
  int ell;
  int cm = *m;
  double cx = *x;
  double coef=sqrt(0.5*M_PI)/pow(-1,cm);
  double pbar1, pbar2;
  //legcoeff has a size of *l
  // the size of pbar, pbar_deri and pbar_deri2 is *l - *m +1
  if(cm>0 && cm<l)
    {
      pbar1 = coef*sqrt((l-cm)*(l+cm+1))*gsl_sf_legendre_sphPlm(l,cm+1,cx);//
      pbar2 = coef*sqrt((l-cm+1)*(l+cm))*gsl_sf_legendre_sphPlm(l,cm-1,cx);
      //      *pbar = coef*(sqrt((l-cm)*(l+cm+1))*pbar1-sqrt((l-cm+1)*(l+cm))*pbar2);
      *pbar = pbar1 - pbar2;
    }
  else if(cm==0)
    {
      if(cm+1>l)
	{
	  pbar1 = 0;
	}
      else
	{
	  pbar1 = coef*sqrt((l-cm)*(l+cm+1))*gsl_sf_legendre_sphPlm(l,cm+1,cx);//
	}
      if(l==0)
	{
	  *pbar = pbar1;
	}
      else
	{
	  *pbar = pbar1+sqrt(0.5*M_PI*(l+1-cm)*(l+cm))*gsl_sf_legendre_sphPlm(l,1-cm,cx);
	}
    }
  else if(cm==l)
    {
      pbar2 = gsl_sf_legendre_sphPlm(l,cm-1,cx);
      *pbar = coef*(-sqrt((l-cm+1)*(l+cm)))*pbar2;
    }
}


void fcal_deri1_array(double * pbar, int * lmax, int * m, double * x)
{
  int ell;
  int cm = *m; int clmax = *lmax;
  for (ell = cm; ell <= clmax; ell++)
    {
      fcal_deri1(&pbar[ell-cm], ell, m, x);
    }   
}
