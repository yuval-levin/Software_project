#include <stdio.h>
#include "spmat.h"

struct graph 
{
  struct spmat_ll A;
  long* vectorDegrees;
  long M;
}