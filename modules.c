#include <stdio.h>
#include "spmat.h"


## RECONSIDER LOCATION OF THOSE STRUCTS, NOTE THEY'RE DEFINED HERE AND IN HEADER - NOT SURE WHERE IS CORRECT
struct graph 
{
  struct spmat_ll A;
  long* vectorDegrees;
  long M;
}

struct division
{
	int len  //number of divisions
	int** divisions //"P" in algorithm
}