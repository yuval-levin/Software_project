#ifndef _POWER_ITER_H
#define _POWER_ITER_H

void create_abk_vec(int rowLength, double* currentB, double* newB,
		struct shiftedDivisionGroup* g, struct graph* graph);

double* create_eigenvector(int rowLength, struct shiftedDivisionGroup* g,
		struct graph* graph);

#endif
