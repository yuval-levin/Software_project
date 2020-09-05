/*
 * Algorithm3.h
 *
 *  Created on: 5 Sep 2020
 *      Author: 97254
 */

#ifndef MODULARITYMAXIMIZATION_H_
#define MODULARITYMAXIMIZATION_H_

double* modularityTimesS(struct graph* graph, int* vectorS,
		struct divisionGroup* g, double sumKiSi);
double dotProduct(double* a, double* b, int col);
void modularityMaximization(struct graph* graph, int* vectorS,
		struct divisionGroup* g);
double sumOfDegreeByVectorS(struct graph* graph, int* vectorS,
		struct divisionGroup* g);

#endif /* MODULARITYMAXIMIZATION_H_ */
