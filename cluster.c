#include<stdio.h>
#include <stdlib.h>
#include "modules.h"

int main(int args, char** argv)
{


}

/* PARAMS: n isn number of nodes */
//SHOULD WE USE TYPEDEF FOR DIVISION?
/* by default: contains 1 divisions, which includes all the nodes.*/
struct division createDivision(int n)
{
	struct division div;

	div=(struct division*)malloc(sizeof(struct division));
	div->len = n;
	div->divisions = (int**)malloc(n*sizeof(int*));
	div[0] = createDefaultDivision(n);

	return div;
}

int* createDefaultDivision(int n)
{
	int i;
	int* vectorOfDivision = (int*)malloc(n*sizeof(int));
	for (i = 0;i < n;i++)
	{
		vectorOfDivision[i] = 1;
	}
	return vectorOfDivision;
}
