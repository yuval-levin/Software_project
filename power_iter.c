#include<stdio.h>
#include <stdlib.h>
#include "modules.h"
#include <math.h>

void createB(double* b, int col) {
	/*creates a random vector */
	int i;
	for (i = 0; i < col; i++) {
		*b = rand();
		b = b + 1;
	}

}

int difference(double * a, double *b, int col) {
	/*returns 1 if: difference of some coordinate in a and b differ by more than epsilon,
	 * 0 otherwise
	 */
	int k;
	double dif = 0;
	double* vec1;
	double* vec2;

	vec1 = a;
	vec2 = b;
	for (k = 0; k < col; k++) {
		dif = fabs(((*vec1) - (*vec2)));
		if (dif >= epsilon) {
			return 1;
		}
		vec1 += 1;
		vec2 += 1;
	}

	return 0;

}

double dotProduct(double* a, double* b, int col) {
	/*dot product of vectors a and b
	 * */
	int k;
	double dot;
	double* vec1;
	double* vec2;
	dot = 0;
	vec1 = a;
	vec2 = b;
	for (k = 0; k < col; k++) {
		dot += ((*vec1) * (*vec2));
		vec1 += 1;
		vec2 += 1;
	}
	return dot;
}

double magnitude(double* vec, int col) {
	/*returns magnitude (norm) of vec with col columns*/
	return sqrt(dotProduct2(vec, vec, col));
}

void divideByMagnitude(double* vec, double magnitude, int col) {
	/*dividing a vector by its' magnitude */
	int i;

	for (i = 0; i < col; i++) {
		*vec = (*vec) / magnitude;
		vec += 1;
	}

}

void updateB(double* b, double* newB, double c) {
	/*given a  vector b, we update its value to be the values of newB */
	int i;
	for (i = 0; i < c; i++) {
		*b = (*newB);
		b++;
		newB++;
	}

}


double spmatProductWithVectorb(int rowIndex, double* vector,
		struct shiftedDivisionGroup* g, struct graph* graph,
		double KjBjDividedByM, double KjDivdedByM) {
	double rowResult = 0;
	double rowSum = 0, ki = graph->vectorDegrees[rowIndex], bi =
			vector[rowIndex];
	struct divisionGroup group = g->group; //eran: consider making it a const; you know no O3 and stuff...

	struct spmat_node* cur_node =get_private(group->groupSubmatrix)[rowIndex];
	rowResult += sumHelper(cur_node, vector, &rowSum); //A times b

	rowResult -= ((KjBjDividedByM) * ki);
	rowResult -= ((rowSum - ki * KjDivdedByM) * bi);
	rowResult += (g->norm * bi);

	return rowResult;

}

//helper function calculates row (given by cur_node) sum , and row times vector v;
// it does two things to prevent iterating the same row twice.
double sumHelper(struct spmat_node* cur_node, double *v, double* rowSum) {
	int index;
	double sum = 0;
	while (cur_node != NULL) {
		index = cur_node->index;
		*rowSum = (*rowSum) + cur_node->data;
		sum += (cur_node->data) * (v[index]);
		cur_node = (struct spmat_node*) cur_node->next;
	}
	return sum;
}

double sparseRowSum(int rowIndex, struct shiftedDivisionGroup* g) {

}

void spmatProductHelperKjBjDividedByM(double* vector,
		struct shiftedDivisionGroup* g, struct graph* graph,
		double* KjBjMultiply, double* KjBj) {
	double sum = 0, sumMultiply = 0;
	int i;
	struct divisionGroup group = g->group;
	for (i = 0; i < group->groupSize; i++) {
		sumMultiply += (vector[i] * graph->vectorDegrees[i]);
		sum += graph->vectorDegrees[i];
	}

	*KjBjMultiply = sumMultiply / (graph->M);
	*KjBj = sum / (graph->M);

}

void createAbkVec( int rowLength, double* currentB, double* newB,
		struct shiftedDivisionGroup* g, struct graph* graph) {
	int i;
	double Abk, KjBjDividedByM, KjDivdedByM;
	int n;

	spmatProductHelperKjBjDividedByM(currentB, g, graph, &KjBjDividedByM,
			&KjDivdedByM); //helper Sums for all rows

	for (i = 0; i < rowLength; i++) {

		/*calculate vector Abk in current coordinate by doing dot prodct of current matrix row with current b vector */
		Abk = spmatProductWithVectorb(i, currentB, g,graph, KjBjDividedByM,
				KjDivdedByM);
		/*updating vector Abk in current coordinate */
		*newB = Abk;
		/*move nextb to next coordinate */
		newB += 1;
	}

}

double* createEigenvalue( int rowLength, struct shiftedDivisionGroup* g,struct graph* graph) {
	//todo, include epsilon for differnece function
	/*since row=col in cov matrix, we use only row param*/
	double* b;
	double* nextb;
	double* covRow;
	double* origNextB;
	int n;
	int dif;

	b = (double*) malloc(rowLength * sizeof(double));
	//TODO: add exit
	createB(b, rowLength);
	covRow = (double*) malloc(rowLength * sizeof(double));
	//TODO: add exit
	nextb = (double*) malloc(rowLength * sizeof(double));
	//TODO: add exit
	/*saving its' original start pointer*/
	origNextB = nextb;
	do {

		createAbkVec( rowLength, b, origNextB, g, graph);
		/*normalizing nextb*/
		divideByMagnitude(origNextB, magnitude(origNextB, rowLength), rowLength);
		/*checking difference between "old" b and "new" b  vectors:*/

		dif = difference(b, origNextB, rowLength);

		/*updating b: */
		updateB(b, origNextB, rowLength);

	} while (dif == 1);

	free(origNextB);
	free(covRow);
	return b;

}


