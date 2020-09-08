#include <stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "ModularityMaximization.h"
#include "OneNorm.h"
#include "power_iter.h"
#include "spmat.h"
#include "Algorithm3.h"


void computeS(double* u1, double* s, int n)
{
	int i;
	for (i = 0; i < n; i++)
	{
		s[i] = u1[i] >= 0 ? 1 : -1;
	}

}


struct shiftedDivisionGroup* newShiftedDivsionGroup(struct divisionGroup* g,struct graph* graph)
{
	struct shiftedDivisionGroup* shiftedG = (struct shiftedDivisionGroup*) malloc(sizeof(struct shiftedDivisionGroup));
	shiftedG->group = g;
	shiftedG->norm = one_norm(graph,g);
	return shiftedG;
}

double computeLeadingEigenvalue(struct shiftedDivisionGroup* shiftedG,double* eigenvector,struct graph* graph)
{
	double eigenvalue,numerator,denominator,rowLength;
	struct divisionGroup* g = shiftedG->group;
	rowLength =g->groupSize;

	double* BShiftedTimesEigenvector = (double*)malloc(rowLength*sizeof(double));
	createAbkVec(rowLength, eigenvector, BShiftedTimesEigenvector,
			shiftedG,graph);
	numerator = dotProduct(BShiftedTimesEigenvector, eigenvector, rowLength);
	denominator = dotProduct(eigenvector,eigenvector,rowLength);
	/*todo: check division by zero and add exit*/
	eigenvalue = numerator / denominator;
	eigenvalue -= shiftedG->norm;

	free(BShiftedTimesEigenvector);
	return eigenvalue;
}

void fillVectorWithOnes(double* vector, int length) {
	int i;
	for (i = 0; i < length; i++) vector[i] = 1;
}

void Algorithm2(double* vectorS, struct divisionGroup* g, struct graph* graph) {
	double eigenvalue, sumKiSi,rightArgument,AtimesS,leftArgument;
	double* eigenvector;

	struct shiftedDivisionGroup* shiftedG = newShiftedDivsionGroup(g,graph);

	eigenvector = createEigenvalue(g->groupSize, shiftedG, graph);
	eigenvalue = computeLeadingEigenvalue(shiftedG,eigenvector,graph);

	if (eigenvalue <= epsilon) {
		/*g is undivisble so S stays the same - everyone are '1';*/
		fillVectorWithOnes(vectorS, g->groupSize);
	}
	else {
		computeS(eigenvector, vectorS, g->groupSize);
		sumKiSi = sumOfDegreeByVectorS(graph, vectorS, g);
		rightArgument = dotProduct(vectorS, modularityTimesS(graph, vectorS, g, sumKiSi), g->groupSize);
		mult_ll(g->groupSubmatrix,vectorS, &AtimesS);
		leftArgument = dotProduct(vectorS,&AtimesS,g->groupSize);
		if (leftArgument+rightArgument<= epsilon)
		{
			/*g is undivisble so S stays the same - everyone are '1';*/
			fillVectorWithOnes(vectorS, g->groupSize);
		}
	}

	free(shiftedG);
}


