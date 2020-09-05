#include<stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "spmat.c"

void flipVectorEntry(int* vector, int entry) {
	vector[entry] = vector[entry] * (-1);
}

/* we wish for S be in the same state it was when we reached maxImproved Score.
 * so we reverse everything that came after it*/
void updateS(int* vectorS, int* indiceVector, int maxImprovedIndex, int length) {
	int i;
	for (i = maxImprovedIndex + 1; i < length; i++) {
		flipVectorEntry(vectorS, indiceVector[i]);
	}
}


void updateImprovedVector(int* improvedVector, int entryIndex, double score) {
	if (entryIndex == 0)
		improvedVector[0] = score;
	else
		improvedVector[entryIndex] = improvedVector[entryIndex - 1] + score;
}

void removeFromUnmoved(struct node* prevOfBiggest, struct node* unmoved) {
	struct node* removedNode;
	if (prevOfBiggest == NULL)
		unmoved = unmoved->next, removedNode = unmoved; //todo: make sure this works
	else
		prevOfBiggest->next = removedNode = prevOfBiggest->next, prevOfBiggest->next->next;
	free(removedNode);
}
double calculateChangeModularity(struct graph* graph, int* vectorS,
		struct divisionGroup* g, double sumKiSi, double prevModularity,
		int changedIndex) {
	double result;
	int degree = graph->vectorDegrees[changedIndex], vectorSChangedIndex =
			vectorS[changedIndex];
	result = prevModularity
			- 4 * vectorSChangedIndex * (degree / graph->M)
					* (sumKiSi - (degree * vectorSChangedIndex));

	return result;
}

//TODO: add explanation.
double* secondArgumentInCalc(struct graph* graph, int* vectorS,
		struct divisionGroup* g, double sumKiSi) {
	int i, index;
	double M = graph->M;
	struct spmat_node* current = g->groupSubmatrix->private[0];
	double* KiDividedByMPlusSum = (double*) malloc(
			g->groupSize * sizeof(double));

	if (KiDividedByMPlusSum == NULL)
		exit(1); //TODO: print error before exit.

	//two iterations are a must, cause we need to find sum first..
	for (i = 0; i < g->groupSize; i++) {
		KiDividedByMPlusSum[i] = (graph->vectorDegrees[current->index] / M)
				* sumKiSi;
		current = current->next;
	}
	return KiDividedByMPlusSum;

}

double* modularityTimesS(struct graph* graph, int* vectorS,
		struct divisionGroup* g, double sumKiSi) {
	int i;
	double* resVec = (double*) malloc(g->groupSize * sizeof(double));
	if (resVec == NULL)
		exit(1); //TODO: print error before exit.
	double* KiDividedByMPlusSum = secondArgumentInCalc(graph, vectorS, g,
			sumKiSi);

	for (i = 0; i < g->groupSize; i++) {
		resVec[i] = (g->sumOfRows[i] * vectorS[i]) + KiDividedByMPlusSum[i];
	}

	free(KiDividedByMPlusSum);

	return resVec;
}


//calculate kisi
double sumOfDegreeByVectorS(struct graph* graph, int* vectorS,
		struct divisionGroup* g) {
	double sum = 0;
	int i;
	struct spmat_node* current = g->groupSubmatrix->private[0];
	for (i = 0; i < g->groupSize; i++) {
		//vectorDegrees is size Of number of nodes in the original A matrix
		sum = sum + (vectorS[i] * graph->vectorDegrees[current->index]); //vectorS is size of g, we use i
		current = current->next;
	}
	return sum;
}

struct node* createUnmovedList(int sizeOfg) {
	int i;
	struct node* head, prev = NULL;
	for (i = 0; i < sizeOfg; i++) {
		prev = appendToList(prev, i);
		if (i == 0)
			head = prev;
	}
	return head;
}

struct node* appendToList(struct node* prev, int index) {
	struct node* current;

	current = (struct node*) malloc(sizeof(struct node*));
	if (current == NULL)
		exit(1); //TODO: print error before exit.
	current->data.num = index;
	current->next = NULL;
	if (prev != NULL)
		prev->next = current;

	return current;
}

//TODO: remove duplicate and add to module of matrix functions
double dotProduct(double* a, double* b, int col) {
	/*dot product of vectors a and b*/
	int k;
	double* vec1 = a, vec2 = b;
	double dot = 0;

	for (k = 0; k < col; k++) {
		dot += ((*vec1) * (*vec2));
		vec1 += 1;
		vec2 += 1;
	}
	return dot;
}


//TODO: is DeltaModularity double int long?
void modularityMaximization(struct graph* graph, int* vectorS,
		struct divisionGroup* g) {

	double modularityChange, Q0, Q1, maxModularityChange, maxImprovedIndex = 0,
			maxImproveScore, sumKiSi;
	int i, indexOfBiggestIncrease, switchFirstUnmovedIteration = 1;
	struct node* unmoved, currentNode, prev, prevOfBiggest;
	double* improvedVector = (double*) malloc(g->groupSize * sizeof(double));
	int* indiceVector = (int*) malloc(g->groupSize * sizeof(int));

	do {
		//improving delta Q by moving ONE index
		unmoved = createUnmovedList(g->groupSize);

		for (i = 0; i < g->groupSize; i++) {
			currentNode = unmoved;
			prev = NULL, prevOfBiggest = NULL;
			sumKiSi = sumOfDegreeByVectorS(graph, vectorS, g);
			if (i == 0)
				Q0 = dotProduct(vectorS,
						modularityTimesS(graph, vectorS, g, sumKiSi));
			else
				Q0 = Q0
						+ calculateChangeModularity(graph, vectorS, g, sumKiSi,
								Q0, indexOfBiggestIncrease);
			//change in Modularity is Q1-Q0 (For Q1 with biggest modularity). we wish to find Q1 so it is Q1-Q0+Q0 ^
			//finding vertex with maximal increase in modularity
			while (currentNode != NULL) {
				flipVectorEntry(vectorS, currentNode->data.num);
				Q1 = calculateChangeModularity(graph, vectorS, g, sumKiSi, Q0,
						currentNode->data.num);

				if (switchFirstUnmovedIteration == 1
						|| Q1 > maxModularityChange) {
					maxModularityChange = Q1, indexOfBiggestIncrease = i;
					maxImproveScore = maxModularityChange;
					prevOfBiggest = prev;
					switchFirstUnmovedIteration = 0;
				}
				flipVectorEntry(vectorS, currentNode->data.num);

				prev = currentNode;
				currentNode = currentNode->next;
			}
			//moving vertex with maximal increase in modularity
			flipVectorEntry(vectorS, indexOfBiggestIncrease);
			removeFromUnmoved(prevOfBiggest, unmoved);
			//updating vector of indice to save the index we now moved:
			indiceVector[i] = indexOfBiggestIncrease;
			//updating the current 'state score'
			updateImprovedVector(improvedVector, i, maxModularityChange); // incrementing scores
			if (improvedVector[i] > maxImproveScore)
				maxImproveScore = improvedVector[i], maxImprovedIndex = i;
			switchFirstUnmovedIteration = 1;

		}
		modularityChange = maxImproveScore;
		updateS(vectorS, indiceVector, maxImprovedIndex, g->groupSize);

	} while (modularityChange > epsilon);

	free(improvedVector);
	free(indiceVector);
	free(unmoved);
}

