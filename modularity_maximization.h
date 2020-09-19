#ifndef _MODULARITY_MAXIMIZATION_H
#define _MODULARITY_MAXIMIZATION_H

double* modularityTimesS(struct graph* graph, double* vectorS,struct divisionGroup* g, double sumKiSi);
double dot_product(double* a, double* b, int col);
void modularityMaximization(struct graph* graph, double* vectorS,struct divisionGroup* g);
double sumOfDegreeByVectorS(struct graph* graph, double* vectorS,struct divisionGroup* g);

#endif
