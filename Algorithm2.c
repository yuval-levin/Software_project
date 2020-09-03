#include<stdio.h>
#include <stdlib.h>
#include "modules.h"

//TODO: is u1 int? long? double?

void computeS(double* u1, int* s, int n) {
	int i;
	for (i = 0; i < n; i++)
		s[i] = u1[i] >= 0 ? 1 : -1;
}

void Algorithm2(int* vectorS, struct divisionGroup* g, struct graph* graph) {
	double eigenvalue, sumKiSi;
	double* eigenvector;

	eigenvector = computeLeadingEigenvector();
	eigenvalue = computeLeadingEigenvalue();
	if (eigenvalue <= epsilon) {
		//g is undivisble so S stays the same - everyone are '1';
		fillVectorWithOnes(vectorS, g->groupSize);
	} else {
		computeS(eigenvector, vectorS, g->groupSize);
		sumKiSi = sumOfDegreeByVectorS(graph, vectorS, g);
		if (dotProduct(vectorS, modularityTimesS(graph, vectorS, g, sumKiSi))
				<= epsilon) {
			//g is undivisble so S stays the same - everyone are '1';
			fillVectorWithOnes(vectorS, g->groupSize);
		}
	}

	//TODO: step 5. Yuval's
}

void fillVectorWithOnes(int* vector, int length) {
	int i;
	for (i = 0; i < length; i++)
		vector[i] = 1;
}

