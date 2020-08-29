#include<stdio.h>
#include <stdlib.h>
#include "modules.h"


/*PARAMS: P,O are the divisions from the algorithm.
		They are created and initialized outside of the function,
		to enable memory freeing easily and outside of the function */
// @CR I think it makes more sense to get the graph and output the division and manage the memory allocations internally.
struct division* Algorithm3(struct division* P, struct division* O)
{

	struct divisionGroup* g; // @CR always give an initial value -probably null in this case. Good practice that will save you lots of debugging.
	short  gWasDivided; //1 for "true", 0 for "false"; // @CR why start with g? What convention are you following? Also in booleans its usually 0 for false and anything else for true.
	while (P->len > 0)
	{
		g = (P->divisions)[0];
		gWasDivided = Algorithm2(&P,0);
		if (gWasDivided == 0) updateOAndP(&P,&O,0); //TODO: if we always call with 0, add 0 inside function instead of argument
		else addToOSize1Groups(&O,&P);
		// STEP 3.4.2 is automatic, since Algorithm2 edits P.
		// So g1,g2 still remain in P is their size is > 0.
	}

	return &O; // @CR bug - should return O;
}
