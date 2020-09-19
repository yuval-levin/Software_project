#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include "modules.h"
#include "spmat.h"
#include "ModularityMaximization.h"
#include "error_codes.h"
#include <time.h>

/*creates a random vector */
void createB(double* b, int col) {

	int i;
	for (i = 0; i < col; i++) {
		*b = rand();
		b = b + 1;
	}
}

/*returns 1 if: difference of some coordinate in a and b differ by more than epsilon,
 * 0 otherwise
 */
int difference(double * a, double *b, int col) {
	int k;
	double dif = 0;
	double* vec1;
	double* vec2;

	vec1 = a;
	vec2 = b;
	for (k = 0; k < col; k++) {
		dif = fabs(((*vec1++) - (*vec2++)));
		if (dif >= epsilon) {
			return 1;
		}
	}
	return 0;
}

/*returns magnitude (norm) of vec with col columns*/
double magnitude(double* vec, int col) {
	return sqrt(dotProduct(vec, vec, col));
}

/*dividing a vector by its' magnitude */
void divideByMagnitude(double* vec, double magnitude, int col) {
	int i;

	for (i = 0; i < col; i++) {
		*vec = (*vec) / magnitude;
		vec += 1;
	}
}

/*given a  vector b, we update its value to be the values of newB */
void updateB(double* b, double* newB, double c) {

	int i;
	for (i = 0; i < c; i++) {
		*b = (*newB);
		b++;
		newB++;
	}
}

/*helper function calculates row (given by cur_node) sum , and row times vector v;
 / it does two things to prevent iterating the same row twice.*/
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

/*helper function
 * Calculates row "rowIndex" of product between sparse matrix (Adjacency matrix A) and b
 * when b is the vector from power_iteration that is constantly changed.
 * */
double spmatProductWithVectorb(int rowIndex, double* vector,
		struct shiftedDivisionGroup* g, struct graph* graph,
		double KjBjDividedByM, double KjDivdedByM, double* AtimesB) {
	double rowResult = 0;
	double rowSum;
	double ki;
	double bi;

	rowSum = 0;
	ki = graph->vectorDegrees[rowIndex];
	bi = vector[rowIndex];

	rowResult += AtimesB[rowIndex]; /*A times b*/
	rowResult -= ((KjBjDividedByM) * ki);
	rowResult -= ((rowSum - ki * KjDivdedByM) * bi);
	rowResult += (g->norm * bi);

	return rowResult;
}

void spmatProductHelperKjBjDividedByM(double* vector,
		struct shiftedDivisionGroup* g, struct graph* graph,
		double* KjBjMultiply, double* KjBj) {
	double sum = 0, sumMultiply = 0;
	double* degreesDividedByM = graph->degreesDividedByM;
	int i;
	struct divisionGroup* group = g->group;
	for (i = 0; i < group->groupSize; i++) {
		sumMultiply += (vector[i] * degreesDividedByM[i]);
		sum += degreesDividedByM[i];
	}
	*KjBjMultiply = sumMultiply;
	*KjBj = sum;
}

/*helper function to product Ab_k,
 * As specified in the power iterator algorithm
 * */
void createAbkVec(int rowLength, double* currentB, double* newB,
		struct shiftedDivisionGroup* g, struct graph* graph) {
	int i;
	double Abk, KjBjDividedByM, KjDivdedByM;
	double* AtimesB;
	spmatProductHelperKjBjDividedByM(currentB, g, graph, &KjBjDividedByM,
			&KjDivdedByM); /*helper Sums for all rows*/
	AtimesB = (double*) malloc(sizeof(double) * g->group->groupSize);
	mult_ll(g->group->groupSubmatrix, currentB, AtimesB);

	for (i = 0; i < rowLength; i++) {

		/*calculate vector Abk in current coordinate by doing dot prodct of current matrix row with current b vector */
		Abk = spmatProductWithVectorb(i, currentB, g, graph, KjBjDividedByM,
				KjDivdedByM, AtimesB);
		/*updating vector Abk in current coordinate */
		*newB = Abk;
		/*move nextb to next coordinate */
		newB += 1;
	}
	free(AtimesB);
}

/*function to calculate eigenvalue for matrix of group g.
 * Calculates eigenvalue using power iteration.
 * */
double* createEigenvalue(int rowLength, struct shiftedDivisionGroup* g,
		struct graph* graph) {
	/*since row=col in cov matrix, we use only row param*/
	double* b;
	double* nextb;
	double* covRow;
	double* origNextB;
	int dif;

	b = (double*) malloc(rowLength * sizeof(double));
	if (b == NULL)
		panic(ERROR_MALLOC_FAILED);

	createB(b, rowLength);
	covRow = (double*) malloc(rowLength * sizeof(double));
	if (covRow == NULL)
		panic(ERROR_MALLOC_FAILED);

	nextb = (double*) malloc(rowLength * sizeof(double));
	if (nextb == NULL)
		panic(ERROR_MALLOC_FAILED);

	/*saving its' original start pointer*/
	origNextB = nextb;
	do {
		createAbkVec(rowLength, b, origNextB, g, graph);
		/*normalizing nextb*/

		divideByMagnitude(origNextB, magnitude(origNextB, rowLength),
				rowLength);
		/*checking difference between "old" b and "new" b  vectors:*/

		dif = difference(b, origNextB, rowLength);

		/*updating b: */
		updateB(b, origNextB, rowLength);

	} while (dif == 1);

	free(origNextB);
	free(covRow);

	return b;
}

