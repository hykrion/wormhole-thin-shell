#ifndef CSV_H_INCLUDED
#define CSV_H_INCLUDED

#include <glib.h>

typedef struct
{
  double x;
  double v;
}TCsv;

void csv_read_data(char *file, GArray** arr);

#endif // CSV_H_INCLUDED
