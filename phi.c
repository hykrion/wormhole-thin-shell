#include "header/phi.h"

#include "header/ui.h"
#include "header/tortoise.h"
#include "header/equation.h"
#include "header/csv.h"

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_deriv.h>
#include <glib.h>

#define SYS_DIM 2
#define PHI_SYS_IND_0 15
#define PHI_SYS_IND_1 20

struct rsSqrtParams
{
  double A;
  double B;
  double rx;
};

struct rxParams
{
  double aa;
};

static double *m_phiRea;
static double *m_phiImg;
static double *m_T;
static double *m_R;
static double *m_v;
static double *m_deltaL;
static int m_nodes;
static GArray* m_csvPhiRea;
static GArray* m_csvPhiIma;
static GArray* m_csvVeff;
/**
  @brief  z(x) para calcular z'(x)
*/
static double
funcZx(double x,
       void *params)
{
  (void)(params);

  /*return sqrt(0.5*(pow(x, 2) + sqrt(pow(x, 2) + 4*pow(rc(), 2))));*/
  return sqrt(0.5*(pow(x, 2) + sqrt(pow(x, 2) + 4)));
}

/**
  @brief  Calcular z'(x)
*/
static double
zx(double x)
{
  gsl_function Fzx;

  double result;
  double absErr;
  double h = 1e-8;

  Fzx.function = &funcZx;
  Fzx.params = NULL;
  gsl_deriv_central (&Fzx, x, h, &result, &absErr);

  return result;
}

static double
funcZyy(double y,
        void *params)
{
  (void)(params);

  double h = 1e-3;
  double plus = A(y + h)*zp(y + h)*zx(y + h);
  double minus = A(y - h)*zp(y - h)*zx(y - h);

  return (plus - minus)/(2*h);
}

/*
  @brief  Ecuación de la barrera/pozo potencial. Es necesario haber calculado
          antes xy[]
*/
static double
phi_v(double x)
{
  (void)x;
  return 0.0;
  /* NOTE No soy capaz de calcularlo...
  double l = ui_get_l();

  /* Para las derivadas *
  gsl_function Fzyy;
  double h = 1e-6;
  double Azpzx;
  double zyy;
  double absErr;
  double xy = tortoise_get_xy_y_indexed(x);

  /* Calcular Azpzx'(y) *
  Fzyy.function = &funcZyy;
  Fzyy.params = NULL;
  gsl_deriv_central(&Fzyy, xy, h, &Azpzx, &absErr);

  zyy = A(xy)*zp(xy)*Azpzx;

  double foo = A(xy);
  double bar = zp(xy);
  double baz = Azpzx;
  double foo1 = (l*(1 + l)*A(xy))/pow(z(xy), 2) + zyy/z(xy);

  /*Veff[y_, l_] := (l (1 + l) A[x]) / z[x]^2 + zyy / z[x] /. x -> TortoiseXY[y];*
  return (l*(1 + l)*A(xy))/pow(z(xy), 2) + zyy/z(xy);
  */
}

/**
  @brief  Sistema de ecuaciones

  @param  r: variable independiente
  @param  y: variables del sistema que permiten tener una ecuación de 1er orden
  @param  dydx: derivadas del sistema: 1a, 2a, 3a...
  @param  params: parámetros para las ecuaciones
*/
static int
funcPhi (double x,
         const double y[],
         double dydx[],
         void* params)
{
  struct TPhiParams *par = (struct TPhiParams *)params;
  double w = (par->w);

  double u = y[0];
  double v = y[1];

  dydx[0] = v;
  dydx[1] = (phi_v(x) - w*w)*u;

  return GSL_SUCCESS;
}

/**
  @brief  Realizar la integración de phi

  @param  IN, a: límite inferior del intervalo
  @param  IN, b: límite superior
  @param  IN, nodes: número de nodos
  @param  IN, ci: condición inicial
  @param  IN, intParams: parámetros necesarios para la integración
  @param  IN, OUT phi: arreglo con los valores de phi calculados

  @return Estado de la integración
*/
static int
phi_integration(double a,
                double b,
                int nodes,
                void *intParams,
                double *phi)
{
  int status = GSL_SUCCESS;

  struct TPhiParams *intPar = (struct TPhiParams*)intParams;
  double w = (intPar->w);
  double u = (intPar->u);
  double v = (intPar->v);

  double h = (b - a)/(nodes - 1);

  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  struct TPhiParams par = {w, 0, 0, 0};

  gsl_odeiv2_system sys = {funcPhi, NULL, SYS_DIM, &par};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd; /* rkf45 rk8pd , rk4imp bsimp msadams msbdf*/
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[2] = {u, v};
  int i;
  for (i = 0; i < nodes; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      printf("x = %f, y = %f\n", x0, y[0]);
      break;
    }
    else
      phi[i] = y[0];
  }
  gsl_odeiv2_driver_free(d);

  return status;
}

/* -----------------------------------------------------------------------------
  PUBLIC
----------------------------------------------------------------------------- */
void phi_init(void)
{
  m_nodes = ui_get_nodes();
  m_phiRea = calloc(m_nodes, sizeof(double));
  m_phiImg = calloc(m_nodes, sizeof(double));
  m_T = calloc(m_nodes, sizeof(double));
  m_R = calloc(m_nodes, sizeof(double));
  m_v = calloc(m_nodes, sizeof(double));
  m_deltaL = calloc(m_nodes, sizeof(double));
  m_csvPhiRea = g_array_sized_new(FALSE, FALSE, sizeof(TCsv), 1024);
  m_csvPhiIma = g_array_sized_new(FALSE, FALSE, sizeof(TCsv), 1024);
  m_csvVeff = g_array_sized_new(FALSE, FALSE, sizeof(TCsv), 1024);
}

void phi_destroy(void)
{
  free(m_phiRea);
  free(m_phiImg);
  free(m_T);
  free(m_R);
  free(m_v);
  free(m_deltaL);
  g_array_free(m_csvPhiRea, FALSE);
  g_array_free(m_csvPhiIma, FALSE);
  g_array_free(m_csvVeff, FALSE);
}

/**
  @brief  Calcula los coeficientes R y T entre un intervalo [a, b] para un rango
          de frecuencias angulares [wMin, wMax] para un momento angular concreto
          ,l, resolviendo la EDO y'(x), que nos servirá una vez resuelta para
          obtener la CI de la EDO x'(y)

  @param  IN, l:  momento angular
*/
void
phi_wave_yx_xy(double l)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double wMin = ui_get_wMin();
  double hW = ui_get_hW();
  int nW = ui_get_nW();
  double w;

  /*
  Puntos alejados de la barrera donde resolveremos el sistema de ecuaciones
  para obtener los valores de R y T
  */
  int i;
  for(i = 0; i < nW; i++)
  {
    /* Calcular 'q' de la ecuación y buscar índices de la barrera */
    w = wMin + hW*i;

    /* Integración hacia atrás */

    /* Parte real */
    double u = cos(w*b);
    double v = -w*sin(w*b);
    struct TPhiParams intParams = {w, l, u, v};
    phi_integration_rea(b, a, m_nodes, &intParams);

    /* Parte imaginaria */
    u = sin(w*b);
    v = w*cos(w*b);
    intParams.u = u;
    intParams.v = v;
    phi_integration_img(b, a, m_nodes, &intParams);

    /* Ordenar los vectores como si hubiéramos integrado hacia delante *
    phi_rea_fwd();
    phi_img_fwd();
    TODO  Peta. Debe ser que el GArray es de 1024, pero de todas formas
          no me hace falta porque ya están ordenados.
    */

    phi_calculate_RT(w, i);
  }
}

/**
  @brief  Realizar la integración de phi

  @param  IN, a: límite inferior del intervalo
  @param  IN, b: límite superior
  @param  IN, nodes: número de nodos
  @param  IN, ci: condición inicial
  @param  IN, intParams: parámetros necesarios para la integración
  @param  IN, OUT phi: arreglo con los valores de phi calculados

  @return Estado de la integración
*/
int
phi_integration_rea(double a,
                    double b,
                    int nodes,
                    void *intParams)
{
  csv_read_data("wave60l0PhiRea.csv", &m_csvPhiRea);
  guint i;
  for(i = 0; i < m_nodes; i++)
  {
    TCsv csv = g_array_index(m_csvPhiRea, TCsv, i);
    m_phiRea[i] = csv.v;
  }
  /*
  NOTE  No soy capaz de calcularlo...
  return phi_integration(a, b, nodes, intParams, m_phiRea)
  */
  return GSL_SUCCESS;
}

/**
  @brief  Realizar la integración de phi

  @param  IN, a: límite inferior del intervalo
  @param  IN, b: límite superior
  @param  IN, nodes: número de nodos
  @param  IN, ci: condición inicial
  @param  IN, intParams: parámetros necesarios para la integración
  @param  IN, OUT phi: arreglo con los valores de phi calculados

  @return Estado de la integración
*/
int
phi_integration_img(double a,
                    double b,
                    int nodes,
                    void *intParams)
{
  csv_read_data("wave60l0PhiIma.csv", &m_csvPhiIma);
  guint i;
  for(i = 0; i < m_nodes; i++)
  {
    TCsv csv = g_array_index(m_csvPhiIma, TCsv, i);
    m_phiImg[i] = csv.v;
  }
  /*
  NOTE  No soy capaz de calcularlo...
  return phi_integration(a, b, nodes, intParams, m_phiImg);
  */
  return GSL_SUCCESS;
}

/**
  @brief  Ordenar phiRea como si se hubiera integrado hacia delante
*/
void
phi_rea_fwd(void)
{
  /* TODO Probar phi_reverse */
  double *tmp = calloc(m_nodes, sizeof(double));

  int i;
  for(i = 0; i < m_nodes; i++)
    tmp[i] = m_phiRea[m_nodes - 1 - i];
  for(i = 0; i < m_nodes; i++)
    m_phiRea[i] = tmp[i];

  free(tmp);
}

/**
  @brief  Ordenar phiRea como si se hubiera integrado hacia delante
*/
void
phi_img_fwd(void)
{
  /* TODO Probar phi_reverse */
  double* tmp = calloc(m_nodes, sizeof(double));

  int i;
  for(i = 0; i < m_nodes; i++)
    tmp[i] = m_phiImg[m_nodes - 1 - i];
  for(i = 0; i < m_nodes; i++)
    m_phiImg[i] = tmp[i];

  free(tmp);
}

/**
  @brief  Calcular los coeficientes R y T de la función de onda
*/
void
phi_calculate_RT(double w,
                 int i)
{
  double a = ui_get_a();
  double h = ui_get_h_forward();
  double B0 = m_phiRea[PHI_SYS_IND_0];
  double B1 = m_phiImg[PHI_SYS_IND_0];
  double B2 = m_phiRea[PHI_SYS_IND_1];
  double B3 = m_phiImg[PHI_SYS_IND_1];

  /* Dejamos las constantes en valor positivo */
  double k1 = -1*(a + PHI_SYS_IND_0*h);
  double k2 = -1*(a + PHI_SYS_IND_1*h);
  gsl_complex z1 = gsl_complex_polar(1, k1*w);
  gsl_complex z2 = gsl_complex_polar(1, k2*w);

  double e0r = GSL_REAL(z1);
  double e0i = GSL_IMAG(z1);
  double e1r = GSL_REAL(z2);
  double e1i = GSL_IMAG(z2);

  double Ar;
  double Ai;
  double Br;
  double Bi;
  gsl_complex AA;
  gsl_complex BB;

  double A[] = { e0r, e0i, e0r, -e0i,
                -e0i, e0r, e0i, e0r,
                 e1r, e1i, e1r, -e1i,
                -e1i, e1r, e1i, e1r };

  double B[] = { B0,
                 B1,
                 B2,
                 B3 };
  gsl_matrix_view m = gsl_matrix_view_array (A, 4, 4);
  gsl_vector_view b = gsl_vector_view_array (B, 4);
  gsl_vector *x = gsl_vector_alloc (4);
  int signum;
  gsl_permutation *permutation = gsl_permutation_alloc (4);
  gsl_linalg_LU_decomp (&m.matrix, permutation, &signum);
  gsl_linalg_LU_solve (&m.matrix, permutation, &b.vector, x);

  Ar = gsl_vector_get(x, 0);
  Ai = gsl_vector_get(x, 1);
  Br = gsl_vector_get(x, 2);
  Bi = gsl_vector_get(x, 3);
  AA = gsl_complex_rect(Ar, Ai);
  BB = gsl_complex_rect(Br, Bi);

  m_T[i] = 1.0 / gsl_complex_abs2(AA);
  m_R[i] = gsl_complex_abs2(BB) / gsl_complex_abs2(AA);

  gsl_permutation_free(permutation);
  gsl_vector_free(x);
}

/**
  @brief  Valores del potencial. Así podemos graficarlo
*/
void
phi_calculate_v(void)
{
  csv_read_data("wave60l0Veff.csv", &m_csvVeff);
  guint i;
  for(i = 0; i < m_nodes; i++)
  {
    TCsv csv = g_array_index(m_csvVeff, TCsv, i);
    m_v[i] = csv.v;
  }

  /* NOTE No soy capaz de calcularlo...
  double a = ui_get_a();
  double h = ui_get_h_forward();

  double x;
  int i;

  for (i = 0; i < m_nodes; i++)
  {
    x = a + i*h;
    m_v[i] = phi_v(x);
  }
  */
}

double*
phi_get_rea(void)
{
  return m_phiRea;
}

double*
phi_get_img(void)
{
  return m_phiImg;
}

double*
phi_get_v(void)
{
  return m_v;
}

double*
phi_get_R(void)
{
  return m_R;
}

double*
phi_get_T(void)
{
  return m_T;
}

/**
  @brief  Calcular delta_l para los diferentes momentos angulares, l.
*/
void
phi_calculate_delta_l(double l)
{
  double wMin = ui_get_wMin();
  double nW = ui_get_nW();
  double hW = ui_get_hW();
  double *R = phi_get_R();
  double w;

  int i;
  for(i = 0; i < nW; i++)
  {
    w = wMin + hW*i;
    m_deltaL[i] += (M_PI*(2*l + 1)*(1 - R[i]))/(w*w);
  }
}

double*
phi_get_deltaL(void)
{
  return m_deltaL;
}
