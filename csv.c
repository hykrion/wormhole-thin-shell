#include "header/csv.h"

#include <stdlib.h>
#include <stdio.h>

static const int CSV_LINES = 1001;

void
csv_read_data(char *file,
              GArray** arr)
{
  FILE *fp = fopen(file, "r");
  char buffer[1024];
  int i = 0;
  while(fgets(buffer, sizeof(buffer), fp) != NULL && i <= CSV_LINES)
  {
    char* token = strtok(buffer,",");
    double x = atof(token);
    token = strtok(NULL, ",");
    double v = atof(token);
    TCsv csv = {x, v};
    g_array_append_val(*arr, csv);
    i++;
  }

  fclose(fp);
}
