#include <stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "spmat.h"

/*
 * calculates column sum of Adjacency matrix in column "column"
 */
double column_sum(struct graph* graph, struct divisionGroup* g, int column) {
	double sum = 0;
	int cnt;
	struct spmat_node* currentNode;

	int* groupMembers = g->groupMembers;
	double* vectorDegrees = graph->vectorDegrees;
	int* sumOfRows = g->sumOfRows;
	double secondArgu = graph->degreesDividedByM[column];
	currentNode = get_private(g->groupSubmatrix)[column];

	/*iterate over all rows*/
	cnt = 0;
	while (currentNode != NULL) { /*since matrix is sparse, we don't know how many rows there are */
		sum = sum + (currentNode->data) - sumOfRows[column];
		sum = sum - (vectorDegrees[groupMembers[cnt++]] * secondArgu);
		currentNode = currentNode->next;
	}

	return sum;
}

/* method to calculate the 1-norm of matrix mat.
 */
double one_norm(struct graph* graph, struct divisionGroup* g) {
	double maxColumn = 0;
	double currentSum = 0;
	int i;
	/* iterate through columns*/
	for (i = 0; i < g->groupSize; i++) {
		currentSum = column_sum(graph, g, i);
		if (maxColumn < currentSum)
			maxColumn = currentSum;
	}

	return maxColumn;
}

