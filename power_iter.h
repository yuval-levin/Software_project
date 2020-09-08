/*
 * power_iter.h
 *
 *  Created on: 5 Sep 2020
 *      Author: 97254
 */

#ifndef POWER_ITER_H_
#define POWER_ITER_H_

void createAbkVec( int rowLength, double* currentB, double* newB,struct shiftedDivisionGroup* g, struct graph* graph);
double* createEigenvalue( int rowLength, struct shiftedDivisionGroup* g,struct graph* graph);


#endif /* POWER_ITER_H_ */
