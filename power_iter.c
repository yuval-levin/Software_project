#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include "modules.h"
#include "error_codes.h"
#include "modularity_maximization.h"
#include "spmat.h"

/*creates a random vector */
static void fill_vector_with_random(double* b, int length) {

	int i;
	for (i = 0; i < length; i++) {
		*b = rand();
		b = b + 1;
	}
}

/*returns 1 if: difference of some coordinate in a and b differ by more than epsilon,
 * 0 otherwise
 */
static int difference_between_vector(double * vec1, double *vec2, int length) {
	int k;
	double dif = 0;

	for (k = 0; k < length; k++) {
		dif = fabs(((*vec1++) - (*vec2++)));
		if (dif >= epsilon) {
			return 1;
		}
	}
	return 0;
}

/*returns magnitude (norm) of vec with col columns*/
static double magnitude(double* vec, int length) {
	return sqrt(dot_product(vec, vec, length));
}

static void divide_by_magnitude(double* vec, double magnitude, int length) {
	int i;

	for (i = 0; i < length; i++) {
		*vec = (*vec) / magnitude;
		vec += 1;
	}
}

/*given a  vector b, we update its value to be the values of newB */
static void update_vector_b(double* b, double* newB, int length) {

	int i;
	for (i = 0; i < length; i++) {
		*b = (*newB);
		b++;
		newB++;
	}
}

/*
 * Calculates row "rowIndex" of product between sparse matrix (Adjacency matrix A) and b
 * when b is the vector from power_iteration that is constantly changed.
 * */
static double spmat_product_with_vector_b(int rowIndex, double* vector,
		struct shiftedDivisionGroup* g, struct graph* graph,
		double KjBjDividedByM, double* AtimesB) {
	double rowResult = 0;
	double rowSum;
	double ki;
	double bi;

	rowSum = g->group->sumOfRows[rowIndex];
	ki = graph->vectorDegrees[g->group->groupMembers[rowIndex]];
	bi = vector[rowIndex];

	rowResult += AtimesB[rowIndex]; /*A times b*/
	rowResult -= ((KjBjDividedByM) * ki);
	rowResult -= ((rowSum) * bi);
	rowResult += (g->norm * bi);

	return rowResult;
}

static void spmat_product_helper_KjBj_divided_by_M(double* vector,
		struct shiftedDivisionGroup* g, struct graph* graph,
		double* KjBjMultiply) {
	double  sumMultiply = 0;
	double* degreesDividedByM = graph->degreesDividedByM;
	int i;
	struct divisionGroup* group = g->group;
	for (i = 0; i < group->groupSize; i++) {
		sumMultiply += (vector[i] * degreesDividedByM[i]);
	}
	*KjBjMultiply = sumMultiply;

}

/*helper function to product Ab_k,
 * As specified in the power iterator algorithm
 * */
void create_abk_vec(int rowLength, double* currentB, double* newB,
		struct shiftedDivisionGroup* g, struct graph* graph) {
	int i;
	double Abk, KjBjDividedByM;
	double* AtimesB;
	spmat_product_helper_KjBj_divided_by_M(currentB, g, graph, &KjBjDividedByM); /*helper Sums for all rows*/
	AtimesB = (double*) malloc(sizeof(double) * g->group->groupSize);
	mult_ll(g->group->groupSubmatrix, currentB, AtimesB);

	for (i = 0; i < rowLength; i++) {

		/*calculate vector Abk in current coordinate by doing dot prodct of current matrix row with current b vector */
		Abk = spmat_product_with_vector_b(i, currentB, g, graph, KjBjDividedByM, AtimesB);
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
double* create_eigenvector(int rowLength, struct shiftedDivisionGroup* g,
		struct graph* graph) {
	/*since row=col in cov matrix, we use only row param*/
	double* b;
	double* nextb;
	double* covRow;
	double* origNextB;
	long loopLimiter;/* to prevent infinite loops*/
	long loopCounter;
	int dif;
	loopCounter = 0;
	loopLimiter = g->group->groupSize * 10000;

	b = (double*) malloc(rowLength * sizeof(double));
	if (b == NULL)
		panic(ERROR_MALLOC_FAILED);

	fill_vector_with_random(b, rowLength);
	covRow = (double*) malloc(rowLength * sizeof(double));
	if (covRow == NULL)
		panic(ERROR_MALLOC_FAILED);

	nextb = (double*) malloc(rowLength * sizeof(double));
	if (nextb == NULL)
		panic(ERROR_MALLOC_FAILED);

	/*saving its' original start pointer*/
	origNextB = nextb;
	do {
		loopCounter++;
		create_abk_vec(rowLength, b, origNextB, g, graph);
		/*normalizing nextb*/

		divide_by_magnitude(origNextB, magnitude(origNextB, rowLength),
				rowLength);
		/*checking difference between "old" b and "new" b  vectors:*/

		dif = difference_between_vector(b, origNextB, rowLength);

		/*updating b: */
		update_vector_b(b, origNextB, rowLength);

		if (loopCounter > loopLimiter)
			panic(ERROR_LOOP_LIMIT_REACHED);
	} while (dif == 1);

	free(origNextB);
	free(covRow);

	return b;
}

