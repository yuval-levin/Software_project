#include<stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "spmat.h"


/*TODO: remove duplicate and add to module of matrix functions*/
double dotProduct(double* a, double* b, int col) {
	/*dot product of vectors a and b*/
	int k;
	double* vec1;
	double* vec2;
	double dot = 0;
	vec1 = a;
	vec2 = b;
	for (k = 0; k < col; k++) {
		dot += ((*vec1) * (*vec2));
		vec1 += 1;
		vec2 += 1;
	}
	return dot;
}


void flipVectorEntry(double* vector, int entry) {
	vector[entry] = vector[entry] * (-1);
}

/* we wish for S be in the same state it was when we reached maxImproved Score.
 * so we reverse everything that came after it*/
void updateS(double* vectorS, int* indiceVector, int maxImprovedIndex, int length) {
	int i;
	for (i = maxImprovedIndex + 1; i < length; i++) {
		flipVectorEntry(vectorS, indiceVector[i]);
	}
}


void updateImprovedVector(double* improvedVector, int entryIndex, double score) {
	if (entryIndex == 0)
		improvedVector[0] = score;
	else
		improvedVector[entryIndex] = improvedVector[entryIndex - 1] + score;
}

struct node* removeFromUnmoved(struct node* prevOfBiggest, struct node* unmoved) {
	/*struct node* removedNode;*/
	if (prevOfBiggest == NULL)
	{
		/*removedNode = unmoved;*/
		unmoved = unmoved->next; /*todo: does this actually change unmoved? */
		/*free(removedNode);  FREE CAUSED MAJOR Bug, handle free of list todo*/
		return unmoved;
	}
	/*else*/
	/*removedNode = prevOfBiggest->next ;*/
		prevOfBiggest->next= prevOfBiggest->next->next;
		/*free(removedNode); FREE CAUSED MAJOR Bug, handle free of list todo*/
		return unmoved;




}
double calculateChangeModularity(struct graph* graph, struct divisionGroup* g, double* vectorS,
	 double sumKiSi, double prevModularity,
		int changedIndex) {
	double result;
	double previousSAS;
	double currentSAS;
	double* prevAS;
	double* currentAS ;
	int size;
	int nodeNum, degree, vectorSChangedIndex;

	size = g->groupSize;

	prevAS = (double*)malloc(size*sizeof(double));
	currentAS = (double*)malloc(size*sizeof(double));
	mult_ll(g->groupSubmatrix,vectorS,prevAS);
	previousSAS = dotProduct(vectorS,prevAS,size);
	/*changing Si to be -Si*/
	flipVectorEntry(vectorS, changedIndex);

	mult_ll(g->groupSubmatrix,vectorS,currentAS);
	currentSAS = dotProduct(vectorS,currentAS,size);
	free(prevAS);
	free(currentAS);

	nodeNum = g->groupMembers[changedIndex];
	degree = graph->vectorDegrees[nodeNum];

	vectorSChangedIndex = vectorS[changedIndex];
	result = prevModularity
			 - 4 * vectorSChangedIndex * (degree / graph->M)
					* (sumKiSi - (degree * vectorSChangedIndex));
	result = result + (currentSAS-previousSAS);
	return result;
}

/*TODO: add explanation.*/
double* secondArgumentInCalc(struct graph* graph,
		struct divisionGroup* g, double sumKiSi) {
	int i;
	double M = graph->M;
	double* KiDividedByMPlusSum;

	KiDividedByMPlusSum = (double*) malloc(
			g->groupSize * sizeof(double));

	if (KiDividedByMPlusSum == NULL)
		exit(1); /*TODO: print error before exit.*/

	/*two iterations are a must, cause we need to find sum first..*/
	for (i = 0; i < g->groupSize; i++) {
		KiDividedByMPlusSum[i] = (graph->vectorDegrees[g->groupMembers[i]] / M)
				* sumKiSi;
	}
	return KiDividedByMPlusSum;

}

double* modularityTimesS(struct graph* graph, double* vectorS,
		struct divisionGroup* g, double sumKiSi) {
	int i;
	double* KiDividedByMPlusSum;
	double* resVec = (double*) malloc(g->groupSize * sizeof(double));
	if (resVec == NULL)
		exit(1); /*TODO: print error before exit.*/
	KiDividedByMPlusSum= secondArgumentInCalc(graph, g,
			sumKiSi);

	for (i = 0; i < g->groupSize; i++) {
		resVec[i] = (g->sumOfRows[i] * vectorS[i]) + KiDividedByMPlusSum[i];
	}

	free(KiDividedByMPlusSum);
	return resVec;
}


/*calculate kisi*/
double sumOfDegreeByVectorS(struct graph* graph, double* vectorS,
		struct divisionGroup* g) {
	double sum = 0;
	int i;
	for (i = 0; i < g->groupSize; i++) { /*we don't know how many are really in group, since it sparse. so we use while */
		/*vectorDegrees is size Of number of nodes in the original A matrix*/
		sum = sum + (vectorS[i] * graph->vectorDegrees[g->groupMembers[i]]); /*vectorS is size of g, we use i*/
	}
	return sum;
}

struct node* appendToList(struct node* prev, int index) {
	struct node* current;

	current = (struct node*) malloc(sizeof(struct node));
	if (current == NULL)
		exit(1); /*TODO: print error before exit.*/
	current->data.num = index;
	current->next = NULL;
	if (prev != NULL)
		prev->next = current;

	return current;
}

struct node* createUnmovedList(int sizeOfg) {
	int i;
	struct node* head = NULL;
	struct node* prev = NULL;
	for (i = 0; i < sizeOfg; i++) {
		prev = appendToList(prev, i);
		if (i == 0)
			head = prev;
	}
	return head;
}



double dotProductInt(int* a, double* b, int col) {
	/*dot product of vectors a and b*/
	int k;
	int* vec1;
	double* vec2;
	double dot = 0;
	vec1 = a;
	vec2 = b;

	for (k = 0; k < col; k++) {
		dot += ((*vec1) * (*vec2));
		vec1 += 1;
		vec2 += 1;
	}
	return dot;
}


/*ODO: is DeltaModularity double int long?*/
void modularityMaximization(struct graph* graph, double* vectorS,
		struct divisionGroup* g) {

	double modularityChange, Q0, Q1, maxModularityChange, maxImprovedIndex = 0,
			maxImproveScore, sumKiSi;
	int i, indexOfBiggestIncrease, switchFirstUnmovedIteration = 1;
	struct node* unmoved;
	struct node* currentNode;
	struct node* prev;
	struct node* prevOfBiggest;
	double* improvedVector;
	int* indiceVector;
	int whileCnt = 0;
	printf("%s","hi");
	improvedVector = (double*) malloc(g->groupSize * sizeof(double));
	 indiceVector = (int*) malloc(g->groupSize * sizeof(int));
	printf("%d",whileCnt++);
	do {
		/*improving delta Q by moving ONE index*/
		unmoved = createUnmovedList(g->groupSize);
		printf("%d",whileCnt);
		whileCnt++;
		for (i = 0; i < g->groupSize; i++) {
			currentNode = unmoved;
			prev = NULL, prevOfBiggest = NULL;
			sumKiSi = sumOfDegreeByVectorS(graph, vectorS, g);
			if (i == 0)
				Q0 = dotProduct(vectorS, modularityTimesS(graph, vectorS, g, sumKiSi),g->groupSize);
			else
				Q0 = Q0

						+ calculateChangeModularity(graph, g,vectorS, sumKiSi,
								Q0, indexOfBiggestIncrease);
			switchFirstUnmovedIteration = 1; /*indicator we need to set i = 0 as currentMax*/
			/*change in Modularity is Q1-Q0 (For Q1 with biggest modularity). we wish to find Q1 so it is Q1-Q0+Q0 ^
			finding vertex with maximal increase in modularity*/
			while (currentNode != NULL) {

				Q1 = calculateChangeModularity(graph,g, vectorS, sumKiSi, Q0,
						currentNode->data.num);

				if (switchFirstUnmovedIteration == 1
						|| Q1 > maxModularityChange) {
					maxModularityChange = Q1;
					indexOfBiggestIncrease = currentNode->data.num; /*finding "k" of biggestIncrease;*/
					maxImproveScore = maxModularityChange;
					prevOfBiggest = prev;
					switchFirstUnmovedIteration = 0;
				}
				flipVectorEntry(vectorS, currentNode->data.num);

				prev = currentNode;
				currentNode = currentNode->next;
			}

			/*moving vertex with maximal increase in modularity*/
			flipVectorEntry(vectorS, indexOfBiggestIncrease);
			unmoved = removeFromUnmoved(prevOfBiggest, unmoved);
			/*updating vector of indice to save the index we now moved:*/
			indiceVector[i] = indexOfBiggestIncrease;
			/*updating the current 'state score'*/
			updateImprovedVector(improvedVector, i, maxModularityChange); /*incrementing scores*/
			if (improvedVector[i] > maxImproveScore)
			{
				maxImproveScore = improvedVector[i];
				maxImprovedIndex = i;
			}

		}
		modularityChange = maxImproveScore;
		updateS(vectorS, indiceVector, maxImprovedIndex, g->groupSize);
		printf("%f",modularityChange)
		printf("%s","x");

	} while (modularityChange > epsilon);
	free(improvedVector);
	free(indiceVector);
	printf("%s \n","done modularity maximization - post free");
	/*free(unmoved);* this was freeing 0x0 so removed. but unmoved is still unfreed ! TODO*/
}

