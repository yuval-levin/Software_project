#include <stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "spmat.h"
/*TODO make sure include c is fine*/

double columnSum(struct graph* graph, struct divisionGroup* g, int column) {
	double sum;
	int i;
	struct spmat_node* currentNode = get_private(g->groupSubmatrix)[0];
	/*iterate over all rows*/
	for (i = 0; i < g->groupSize; i++) {
		sum = sum + (currentNode->data) - g->sumOfRows[column];
		sum = sum
				- ((graph->vectorDegrees[i] * graph->vectorDegrees[column])
						/ graph->M);
		currentNode = currentNode->next;
	}
	return sum;
}


/* method to calculate the 1-norm of matrix mat.
 TODO: check if double is necessary */
double one_norm(struct graph* graph, struct divisionGroup* g) {
	double maxColumn = 0, currentSum;
	int i;
	/* iterate through columns*/
	for (i = 0; i < g->groupSize; i++) {
		currentSum = columnSum(graph, g, i);
		if (maxColumn < currentSum)
			maxColumn = currentSum;
	}
	return maxColumn;
}



