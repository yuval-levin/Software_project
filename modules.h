#ifndef HEADER_FILE
#define HEADER_FILE

struct graph
{
  struct spmat_ll A;
  long* vectorDegrees;
  long M;
};

struct division
{
	int len;  //number of divisions
	int** divisions; //"P" in algorithm
};

#endif
