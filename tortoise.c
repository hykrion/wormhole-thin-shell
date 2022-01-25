#include "header/tortoise.h"

#include "header/ui.h"
#include "header/equation.h"
#include "header/csv.h"

#include <math.h>
#include <glib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#define TOR_SYS_DIM 1
#define TOR_MIN_INF -1000000

static GArray* m_csv;
static double *m_xy;
static double *m_yx;
static double *m_xyAnalytical;
static int m_nodes;

/**
  @brief  ODE IV de 1er grado para resolver x[y]' que aparece en el cálculo del
          potencial V.
  @param  t     Variable independiente
  @param  y[]   Parte izq del sistema de ecuaciones de primer grado
  @param  sys[] Parte dcha del sistema de ecuaciones de primer grado
*/
static int
funcXY(double x,
       const double y[],
       double sys[],
       void* params)
{
  (void)(x);
  (void)(params);
  double xy = y[0];

  sys[0] = A(xy)*zp(xy);

  return GSL_SUCCESS;
}

/* -----------------------------------------------------------------------------
  PUBLIC
----------------------------------------------------------------------------- */
void
tortoise_init(void)
{
  m_nodes = ui_get_nodes();
  m_xy = calloc(m_nodes, sizeof(double));
  m_yx = calloc(m_nodes, sizeof(double));
  m_xyAnalytical = calloc(m_nodes, sizeof(double));
  m_csv = g_array_sized_new(FALSE, FALSE, sizeof(TCsv), 1024);
}

void tortoise_destroy(void)
{
  free(m_xy);
  free(m_yx);
  free(m_xyAnalytical);
  g_array_free(m_csv, FALSE);
}

/**
  @brief  Calcular las coordenadas 'yx' y 'xy'
*/
void
tortoise_calculate_yx_xy(void)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double rS = ui_get_rS();
  double aa = 0.0;
  struct tortoise_xyParams torParam = {rS, aa};

  double ci = b - rS*log(b) + pow(rS, 2)/b - 1/pow(b, 2);
  /* Sabemos que en +inf r = r* */
  ci = b;
  tortoise_xy_integration(b, a, m_nodes, ci, &torParam);
  /* NOTE Como no la calculamos sino que viene de Wolfram, no hace falta revertir
  tortoise_reverse_xy();
  */
}

int
tortoise_xy_integration(double a,
                        double b,
                        int nodes,
                        double ci,
                        void *param)
{
  int status = GSL_SUCCESS;

  int n = nodes - 1;

  csv_read_data("t60.csv", &m_csv);
  guint i;
  for(i = 0; i < m_csv->len; i++)
  {
    TCsv csv = g_array_index(m_csv, TCsv, i);
    m_xy[i] = csv.v;
  }
  /*
  NOTE  No consigo integrarlo. Sobre el punto 605 y[0] es -inf
  Como podemos integrar hacia delante o hacia atrás dependiendo de 'a' y 'b',
  calculamos 'h'
  *
  double h = (b - a)/ n;
  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  gsl_odeiv2_system sys = {funcXY, NULL, TOR_SYS_DIM, NULL};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd; /* rk2 rk4 rkf45 rk8pd, msbdf msadams rk4imp *
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[1] = {ci};
  int i;
  for (i = 0; i < nodes; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      break;
    }
    /*
    y[0] valor de la función
    y[1] valor de la derivada de la función (en este caso no existe)
    *
    m_xy[i] = y[0];
  }
  gsl_odeiv2_driver_free(d);
  */
  return status;
}

/**
  @brief  Invertir el arreglo. Útil tras integrar hacia atrás

  TODO    Seguro que hay una forma eficiente de hacer esto
*/
void
tortoise_reverse_xy(void)
{
  double *tmp = calloc(m_nodes, sizeof(double));

  int i;
  for(i = 0; i < m_nodes; i++)
    tmp[i] = m_xy[m_nodes - 1 - i];
  for(i = 0; i < m_nodes; i++)
    m_xy[i] = tmp[i];

  free(tmp);
}

double*
tortoise_get_xy(void)
{
  return m_xy;
}

double*
tortoise_get_yx(void)
{
  return m_yx;
}

double
tortoise_get_xy_i(int i)
{
  return m_xy[i];
}

double
tortoise_get_yx_i(int i)
{
  return m_yx[i];
}

/**
  @brief  Buscar el valor 'x' en 'yx' de un vector ordenado e indizado
*/
double
tortoise_get_xy_y_indexed(double x)
{
  double foo;
  double a = ui_get_a();
  double b = ui_get_b();
  double h = ui_get_h_forward();

  if(x < a)
    x = a;
  else if(x > b)
    x = b;

  int i = (x - a) / h;

  return m_xy[i];
}
