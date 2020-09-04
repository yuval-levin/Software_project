#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

void createB(double* b, int col) {
	/*creates a random vector */
	int i;
	for (i = 0; i < col; i++) {
		*b = rand();
		b = b + 1;
	}

}

int difference(double * a, double *b, int col, double epsilon) {
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

void createAbkVec(double* covRow, int row, double* currentB, double* newB,
		struct shiftedDivisionGroup* g,struct graph* graph) {
	int i;
	double Abk, KjBjDividedByM,KjDivdedByM;
	int n;

	spmatProductHelperKjBjDividedByM(currentB,g,graph,&KjBjDividedByM,&KjDivdedByM);
	for (i = 0; i < row; i++) {
		/* we read the cov matrix one row at a time:*/
		covRow = g->groupSubmatrix.getRow(i)
		n = fread(covRow, sizeof(double), row, input);
		assert(n == row);

		/*calculate vector Abk in current coordinate by doing dot prodct of current matrix row with current b vector */
		Abk = spmatProductWithVectorb(covRow, currentB, row);
		/*updating vector Abk in current coordinate */
		*newB = Abk;
		/*move nextb to next coordinate */
		newB += 1;
	}

}
double spmatProductWithVectorb(int rowIndex, double* vector,
		struct shiftedDivisionGroup* g, struct graph* graph) {

	struct divisionGroup group = g->group; //eran: consider making it a const; you know no O3 and stuff...

	struct spmat_node* cur_node = group->groupSubmatrix->private[rowIndex];
	double Atimesb = sumHelper(cur_node, vector);

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
double* createEigenvalue(FILE* input, int row, double epsilon) {
	/*since row=col in cov matrix, we use only row param*/
	double* b;
	double* nextb;
	double* covRow;
	double* origNextB;
	int n;
	int dif;

	b = (double*) malloc(row * sizeof(double));
	//TODO: add exit
	createB(b, row);
	covRow = (double*) malloc(row * sizeof(double));
	//TODO: add exit
	nextb = (double*) malloc(row * sizeof(double));
	//TODO: add exit
	/*saving its' original start pointer*/
	origNextB = nextb;
	do {

		createAbkVec(covRow, row, b, origNextB, input);
		/*normalizing nextb*/
		divideByMagnitude(origNextB, magnitude(origNextB, row), row);
		/*checking difference between "old" b and "new" b  vectors:*/

		dif = difference(b, origNextB, row, epsilon);

		/*updating b: */
		updateB(b, origNextB, row);

	} while (dif == 1);

	free(origNextB);
	free(covRow);
	return b;

}

int main(int args, char** argv) {
	FILE* input;
	FILE* output;
	double epsilon;
	double* eigenVector;
	int row, n, x;

	//TODO: include modules for epsilon
	/*writing dimnesion of vector: */
	/*constant 1*/

	eigenVector = createEigenvalue(input, row, epsilon);

	free(eigenVector);

	return 0;
}
/*
 * eigen.c
 *
 *  Created on: 13 Apr 2020
 *      Author: 97254
 */

