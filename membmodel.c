#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "c99.h"
#include "name.h"
//#include "fail.h"
//#include "types.h"
//#include "mem.h"
//#include "tensor.h"
//#include "poly.h"
//#include "lob_bnd.h"

# define fmatlaw FORTRAN_NAME(matlaw,MATLAW)

extern struct{
	double hf_gscell;
} fmatlaw;

#if defined(NONHOOK)
# define cal_pwspi_nh FORTRAN_NAME(calpwspi,CALPWSPI)
# define cal_pwspi_deri_nh FORTRAN_NAME(calpwspi_deri,CALPWSPI_DERI)
# define cal_energy_nh FORTRAN_NAME(calenergy,CALENERGY)
#elif defined(SKALAK)
# define cal_pwspi_sk FORTRAN_NAME(calpwspi,CALPWSPI)
# define cal_pwspi_deri_sk FORTRAN_NAME(calpwspi_deri,CALPWSPI_DERI)
# define cal_energy_sk FORTRAN_NAME(calenergy,CALENERGY)
#endif

void cal_pwspi_nh(double * pWspI, double * GsCELL, double * CSK, double * I1, double * I2 )
{
  pWspI[0] =  fmatlaw.hf_gscell;
  pWspI[1] = -fmatlaw.hf_gscell*pow((*I2)+1,-2);
  //  printf("%f %f %f %f %f \n",*(pWspI+1),*GsCELL,*CSK,*I1,*I2);
}

void cal_pwspi_sk(double * pWspI, double * GsCELL, double * CSK, double * I1, double * I2 )
{
  pWspI[0] = fmatlaw.hf_gscell*(1+(*I1));
  pWspI[1] = fmatlaw.hf_gscell*((*CSK) * (*I2) - 1);
}


void cal_energy_nh(double * eshear, double * edilation, double * GsCELL, double * CSK, double * I1, double * I2 )
{
  *eshear    =  fmatlaw.hf_gscell*(*I1-1+1.0/(*I2+1));
  *edilation = 0;
}

void cal_energy_sk(double * eshear, double * edilation, double * GsCELL, double * CSK, double * I1, double * I2 )
{
  double I2sq = (*I2)*(*I2);
  *eshear    =  0.5*fmatlaw.hf_gscell*((*I1)*(*I1) + 2*(*I1) - 2*(*I2) - 0.5*I2sq);
  *edilation = 0.25*(2*(*CSK)+1)*fmatlaw.hf_gscell*I2sq;
}




void cal_pwspi_deri_nh(double * pWspI1,      double * pWspI2,
                       double * pWspI1_deri, double * pWspI2_deri,
                       double * GsCELL, double * CSK, 
                       double * I1, double * I2,
                       double * I1de, double * I2de)
{
  //pWspI1_deri is the derivative of pWspI1 with respect to \xi_{1} and \xi_{2}.
  //I1de is the derivative of I1 with respect to \xi_{1} and \xi_{2}.
  double coeff;
  int al;
  *pWspI1 = fmatlaw.hf_gscell;
  *pWspI2 = -fmatlaw.hf_gscell*pow((*I2)+1,-2);
  coeff   = (*GsCELL)*pow((*I2)+1,-3);
  //  printf("%f \n"),*pWspI1;
  for (al=0;al<2;al++)
    {
      pWspI1_deri[al] = 0;
      pWspI2_deri[al] = coeff*I2de[al];
    }
}

void cal_pwspi_deri_sk(double * pWspI1,      double * pWspI2,
                       double * pWspI1_deri, double * pWspI2_deri,
                       double * GsCELL, double * CSK, 
                       double * I1, double * I2,
                       double * I1de, double * I2de)
{
  //pWspI1_deri is the derivative of pWspI1 with respect to \xi_{1} and \xi_{2}.
  double coeff1,coeff2;
  int al;
  *pWspI1 = fmatlaw.hf_gscell*(1+(*I1));
  *pWspI2 = fmatlaw.hf_gscell*((*CSK) * (*I2) - 1);

  coeff1   = fmatlaw.hf_gscell;
  coeff2   = (*CSK)*coeff1;
  for (al=0;al<2;al++)
    {
      pWspI1_deri[al] = coeff1 * I1de[al];
      pWspI2_deri[al] = coeff2 * I2de[al];
    }
}
