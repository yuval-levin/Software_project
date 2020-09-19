#include <stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "ModularityMaximization.h"
#include "OneNorm.h"
#include "power_iter.h"
#include "spmat.h"
#include "Algorithm3.h"
#include "error_codes.h"

/*helper function:
 * given an vector u1 (which will be our eigenvector)
 * We calculate the S vector created by it,
 * when a positive value (entry of u1 > epsilon) means S's entry will be 1, -1 otheriwse.
 */
static void compute_s(double* u1, double* s, int n) {
	int i;
	for (i = 0; i < n; i++) {
		s[i] = u1[i] >= epsilon ? 1 : -1; /*bigger than 0 is bigger than epsilon*/
	}
}

/*helper function for power iteration process:
 *given a divisionGroup g, we calculate the 1-norm of the Adjacency matrix and create a new divisionGroup,
 *Shifted by the 1-norm calculated - shiftedDivisionGroup
 * */
static struct shiftedDivisionGroup* new_shifted_divsion_group(struct divisionGroup* g,
		struct graph* graph) {
	struct shiftedDivisionGroup* shiftedG =
			(struct shiftedDivisionGroup*) malloc(
					sizeof(struct shiftedDivisionGroup));
	if (shiftedG == NULL)
		panic(ERROR_MALLOC_FAILED);
	shiftedG->group = g;
	shiftedG->norm = one_norm(graph, g);
	return shiftedG;
}

/*helper function for calculation eigenvalue
 * given the eigenvector and the shifted matrix
 * */
static double eigen_value_calc_helper(struct shiftedDivisionGroup* shiftedG,
		double* eigenvector, int rowLength, double* BShiftedTimesEigenvector) {
	double numerator, denominator, eigenvalue;

	numerator = dotProduct(BShiftedTimesEigenvector, eigenvector, rowLength);
	denominator = dotProduct(eigenvector, eigenvector, rowLength);

	if (denominator < epsilon)
		panic(ERROR_DIVISION_BY_ZERO);

	eigenvalue = numerator / denominator;
	eigenvalue -= shiftedG->norm;

	return eigenvalue;
}

/*given a shiftedMatrix, calculates its eigenValue
 * using power iteration
 * */
static double computeLeadingEigenvalue(struct shiftedDivisionGroup* shiftedG,
		double* eigenvector, struct graph* graph) {
	double eigenvalue, rowLength;
	struct divisionGroup* g = shiftedG->group;
	double* BShiftedTimesEigenvector;
	rowLength = g->groupSize;

	BShiftedTimesEigenvector = (double*) malloc(rowLength * sizeof(double));
	if (BShiftedTimesEigenvector == NULL)
		panic(ERROR_MALLOC_FAILED);

	create_abk_vec(rowLength, eigenvector, BShiftedTimesEigenvector, shiftedG,
			graph);

	eigenvalue = eigen_value_calc_helper(shiftedG, eigenvector, rowLength,
			BShiftedTimesEigenvector);
	free(BShiftedTimesEigenvector);
	return eigenvalue;
}

/*helper function, fills given vector with 1.
 * This helps when in Algorithm2, group is undivisble - meaning
 * Vector S should be all "1"
 * */
static void fill_vector_with_ones(double* vector, int length) {
	int i;
	for (i = 0; i < length; i++)
		vector[i] = 1;
}

/*helper function to calculate change in modularity,
 * As specified in Algorithm2 row 4.
 *
 */
static double calc_mod_change(double* vectorS, struct divisionGroup* g,
		struct graph* graph) {
	double sumKiSi;
	double rightArgument;
	double leftArgument;
	double* rightArguTemp;
	double* AtimesS;

	AtimesS = (double*) malloc(g->groupSize * (sizeof(double)));
	if (AtimesS == NULL)
		panic(ERROR_MALLOC_FAILED);

	sumKiSi = sumOfDegreeByVectorS(graph, vectorS, g);
	rightArguTemp = modularityTimesS(graph, vectorS, g, sumKiSi);
	rightArgument = dotProduct(vectorS, rightArguTemp, g->groupSize);
	mult_ll(g->groupSubmatrix, vectorS, AtimesS);
	leftArgument = dotProduct(vectorS, AtimesS, g->groupSize);

	free(rightArguTemp);
	free(AtimesS);

	return leftArgument + rightArgument;
}

/*
 * Algorithm2,
 * As described in sp_project.pdf
 *  */
void Algorithm2(double* vectorS, struct divisionGroup* g, struct graph* graph) {
	double eigenvalue = 0;
	double* eigenvector;
	double changeInModularity;
	struct shiftedDivisionGroup* shiftedG;

	shiftedG = new_shifted_divsion_group(g, graph);

	eigenvector = create_eigenvector(g->groupSize, shiftedG, graph);
	eigenvalue = computeLeadingEigenvalue(shiftedG, eigenvector, graph);

	if (eigenvalue <= epsilon) {
		/*g is undivisble so S stays the same - everyone are '1';*/
		fill_vector_with_ones(vectorS, g->groupSize);
	} else {

		compute_s(eigenvector, vectorS, g->groupSize);
		changeInModularity = calc_mod_change(vectorS, g, graph);

		if (changeInModularity <= epsilon) {
			/*g is undivisble so S stays the same - everyone are '1';*/
			fill_vector_with_ones(vectorS, g->groupSize);
		}

	}
	free(eigenvector);
	free(shiftedG);
}

