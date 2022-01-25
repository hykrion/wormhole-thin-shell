#include "header/equation.h"

#include "header/ui.h"
#include "header/csv.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl_sf.h>
#include <gsl_odeiv2.h>

#define EQU_SYS_DIM 1

static GArray* m_csv;
/* -----------------------------------------------------------------------------
  PRIVATE
----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
  PUBLIC
----------------------------------------------------------------------------- */
void equation_init(void)
{
  /*
  NOTE  Por lo visto lo mejor es inicializar en potencias de 2, y como nosotros
        tenemos 1001 líneas, usaremos 1024
  */
  m_csv = g_array_sized_new(FALSE, FALSE, sizeof(TCsv), 1024);
}

void equation_destroy(void)
{
  g_array_free(m_csv, FALSE);
}

/**
  En realiad es una constante... d_c = 0.5720698226263603;
*/
double
d_c(void)
{
  /*
  double nBI = pow(M_PI, 1.5)/(3*pow(gsl_sf_gamma(0.75), 2));

  return 1.0/(sqrt(2)*nBI);
  */
  return 0.5720698226263603;
}

/**
  NOTE  Me desconcierta...
*/
double
rc(void)
{
  return 1.0;
}

/**
  @brief  r^2(x) = 1/2*(x^2 + sqrt(x^4 + 4*rc^4)) Eq.6
*/
double
r(double x)
{
  /*return sqrt(0.5*(pow(x, 2) + sqrt(pow(x, 2) + 4*pow(rc(), 2))));*/
  return sqrt(0.5*(pow(x, 2) + sqrt(pow(x, 2) + 4)));
}

/**
  @brief  z(x) = r(x)/rc  Eq. 5
*/
double
z(double x)
{
  /*return r(x)/rc();*/
  return r(x);
}

/**
  @brief  z+(x) = 1 + 1/z^4(x)  Eq. 5
*/
double
zp(double x)
{
  return 1.0 + 1.0/pow(z(x), 4);
}

/**
  @brief  z-(x) = 1 - 1/z^4(x)  Eq. 5
*/
double
zm(double x)
{
  return 1.0 - 1.0/pow(z(x), 4);
}

/**
  NOTE  Me desconcierta...
*/
double
rSbyrc(void)
{
  double ns = ui_get_ns();

  return ns/d_c();
}

/**
  @brief  H(x) = -1/d_c + 1/2*sqrt(z^4(x) - 1)*(f_lamb(3/4)(x) + f_lamb(7/4)(x))
          Eq. 8
          Simplificación en Eq. 36-45 paper Reissern-Nordström Black Holes
          in Extended Palatini Theories (1207.6004)
  NOTE    Como no soy capaz de realizar bien el cálculo, voy a usar Mathematica
          para hacer los cálculos en el rango [-60, 60]
*/
double
H(double x)
{
  double result;

  static gboolean fileReaded = FALSE;

  if(!fileReaded)
  {
    fileReaded = TRUE;
    csv_read_data("h60.csv", &m_csv);
  }
  double a = ui_get_a();
  double h = ui_get_h_forward();
  int i = (x - a) / h;

  TCsv csv = g_array_index(m_csv, TCsv, i);
  result = csv.v;

  return result;
}

/**
  @brief  A(x) = 1/z+*(1 - rS/rc*(1 + d_1*H(x))/(z(x)*z-^1/2)) Eq.4
*/
double
A(double x)
{
  /* A[x_] = (1 - rSbyrc (1 + dc H[x]) / (z Sqrt[Zm[x]])) / Zp[z] /. z -> z[x]; */
  return (1.0 - rSbyrc()*(1.0 + d_c()*H(x)) / (z(x)*sqrt(zm(x)))) / zp(z(x));
}
