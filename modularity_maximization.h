#ifndef _MODULARITY_MAXIMIZATION_H
#define _MODULARITY_MAXIMIZATION_H

double* modularity_times_s(struct graph* graph, double* vectorS,struct divisionGroup* g, double sumKiSi);
double dot_product(double* a, double* b, int col);
void modularity_maximization(struct graph* graph, double* vectorS,struct divisionGroup* g);
double sum_of_degree_by_vector_s(struct graph* graph, double* vectorS,struct divisionGroup* g);

#endif
