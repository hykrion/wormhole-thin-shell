#ifndef EQUATION_H_INCLUDED
#define EQUATION_H_INCLUDED

void equation_init(void);
void equation_destroy(void);

double equation_gz_integration(void);
double equation_gz_aprox(void);

double d_c(void);
double rc(void);
double r(double x);
double z(double x);
double zp(double x);
double zm(double x);
double rSbyrc(void);
double H(double x);
double A(double x);


#endif // EQUATION_H_INCLUDED
