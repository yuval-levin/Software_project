#include<stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "spmat.h"
#include "error_codes.h"
#include <time.h> /*todo: remove time etc*/


double calcAiSi (double* vectorS, int changedIndex, struct _spmat* mat){
	struct spmat_node* cur_node;
	cur_node = get_private(mat)[changedIndex];

	return sumTimesVectorHelper(cur_node, vectorS);
}

/*todo: delete this*/
void print_array(double *vec, int dim) {
	int i;
	for (i = 0; i < dim; i++) {
		setbuf(stdout, 0);
		printf("%f ", vec[i]);
	}
	printf("\n");
}
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

/* Helper function:
 * We wish for S be in the same state it was when we reached maxImproved Score.
 * so we reverse everything that came after it*/
void updateS(double* vectorS, int* indiceVector, int maxImprovedIndex,
		int length) {
	int i;
	for (i = maxImprovedIndex + 1; i < length; i++) {
		flipVectorEntry(vectorS, indiceVector[i]);
	}
}

/* As seen in psuedo-code for Algo4 row 14,
 * We update "improvedVector" by Aggregating previous improvements with current improvement.
 * We also keep track of index of maximum improvement !
 * */
void updateImprovedVector(double* improvedVector, int entryIndex, double score,
		double* maxImprovedIndex, double* curMax) {
	/*first improvement is always maximal */
	if (entryIndex == 0) {
		improvedVector[0] = score;
		*curMax = score;
		*maxImprovedIndex = 0;
	} else {
		improvedVector[entryIndex] = improvedVector[entryIndex - 1] + score;
		/*if current improvement is bigger than maximal improvement, we update the max value and index.N*/
		if (improvedVector[entryIndex] > *curMax) {
			*curMax = improvedVector[entryIndex];
			*maxImprovedIndex = entryIndex;
		}
	}

}

struct node* addToList(struct node* list, struct node* node) {
	if (list == NULL) {
		list = node;
		list->next = NULL;
	} else {
		node->next = list;
		list = node;
	}
	return list;
}
/*
 * Helper function to remove node from LinkedList representing "UNMOVED" nodes in algo4.
 * (we remove nodes that had the biggest change in modularity)
 * Since LinkedList is one-directional, removal of node is done by using the previous node to it,
 * and setting its value to null.
 * If previousNode = NULL, then we are removing head of LinkedList.
 * We return the new LinkedList, sans the remove node.
 * */
struct node* removeFromUnmoved(struct node* prevOfBiggest, struct node* unmoved,
		struct node** removedNode) {

	if (prevOfBiggest == NULL) {
		*removedNode = unmoved;
		unmoved = unmoved->next;
	} else {
		*removedNode = prevOfBiggest->next;
		prevOfBiggest->next = prevOfBiggest->next->next;
	}

	return unmoved;
}


/*
 * Helper function:
 * calculates vectorS times A(Adjacency matrix) times vectorS
 * */
double calcSAS(struct divisionGroup* g, double* vectorS) {

	int size;
	double* prevAS;
	double previousSAS;

	size = g->groupSize;
	prevAS = (double*) malloc(size * sizeof(double));

	if (prevAS == NULL)
		panic(ERROR_MALLOC_FAILED);
	/*	flipVectorEntry(vectorS, changedIndex);*/
	mult_ll(g->groupSubmatrix, vectorS, prevAS);
	previousSAS = dotProduct(vectorS, prevAS, size); /*calc SAS prev*/
	/*flipVectorEntry(vectorS, changedIndex);*/
	free(prevAS);

	/*printf("calculate SAS took %f seconds\n", ((double) (end - start) / CLOCKS_PER_SEC));*/

	return previousSAS;
}

/*Calculates Change in Modularity using previous SAS*/
double calculateChangeModularityWithPrevSas(struct graph* graph,
		struct divisionGroup* g, double* vectorS, double sumKiSi,
		double prevModularity, int changedIndex, double previousSAS) {
	double newModularityY;
	double currentSAS;
	double* currentAS;
	int size;
	int nodeNum, degree, vectorSChangedIndex;
	/* vectorS arrives FLIPPED */
	size = g->groupSize;

	currentAS = (double*) malloc(size * sizeof(double));
	if (currentAS == NULL)
		panic(ERROR_MALLOC_FAILED);

	mult_ll(g->groupSubmatrix, vectorS, currentAS);
	currentSAS = dotProduct(vectorS, currentAS, size); /* calc SAS current */

	nodeNum = g->groupMembers[changedIndex];
	degree = graph->vectorDegrees[nodeNum];

	vectorSChangedIndex = vectorS[changedIndex]; /* entry value AFTER FLIP */

	newModularityY = prevModularity
			- 4 * vectorSChangedIndex * (graph->degreesDividedByM[nodeNum])
					* (sumKiSi + (degree * vectorSChangedIndex));
	newModularityY = newModularityY + (currentSAS - previousSAS);

	free(currentAS);

	return newModularityY - prevModularity;
}

double calculateChangeModularity(struct graph* graph, struct divisionGroup* g,
		double* vectorS, double sumKiSi, double prevModularity,
		int changedIndex, double* prevSAS) {
	clock_t start, end;

	double previousSAS, sumAiSi, newModularityY, currentSAS;
	int nodeNum, degree, vectorSChangedIndex;

	start = clock();
	srand(time(NULL));
	previousSAS = *prevSAS;

	nodeNum = g->groupMembers[changedIndex];
	degree = graph->vectorDegrees[nodeNum];

	vectorSChangedIndex = vectorS[changedIndex]; /* entry value AFTER FLIP*/

	/* calc new SAS*/
	sumAiSi = calcAiSi(vectorS, changedIndex, g->groupSubmatrix);			/*TODO: delete*/
	currentSAS = previousSAS - 4 * vectorSChangedIndex * sumAiSi;

	*prevSAS = currentSAS; /*update SAS*/

	newModularityY = prevModularity
			- 4 * vectorSChangedIndex * (degree / graph->M)
					* (sumKiSi + (degree * vectorSChangedIndex));
	newModularityY = newModularityY + (currentSAS - previousSAS);

	return newModularityY - prevModularity;

	end = clock();
	printf("calculate change modularity took %f seconds\n",
			((double) (end - start) / CLOCKS_PER_SEC));

}

/*TODO: add explanation.*/
double* secondArgumentInCalc(struct graph* graph, struct divisionGroup* g,
		double sumKiSi) {
	int i;
	double* KiDividedByMPlusSum;

	KiDividedByMPlusSum = (double*) malloc(g->groupSize * sizeof(double));

	if (KiDividedByMPlusSum == NULL)
		panic(ERROR_MALLOC_FAILED);
	/*two iterations are a must, cause we need to find sum first..*/
	for (i = 0; i < g->groupSize; i++) {
		KiDividedByMPlusSum[i] = (graph->degreesDividedByM[i])
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
		panic(ERROR_MALLOC_FAILED);

	KiDividedByMPlusSum = secondArgumentInCalc(graph, g, sumKiSi);

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

struct node* appendToInitUnmoved(struct node* prev, int index) {
	struct node* current;

	current = (struct node*) malloc(sizeof(struct node));
	if (current == NULL)
		panic(ERROR_MALLOC_FAILED);

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
		prev = appendToInitUnmoved(prev, i);
		if (i == 0)
			head = prev;
	}
	return head;
}

/* remove  todo*/
void printG(struct divisionGroup* g) {
	int n = g->groupSize;
	int i;
	struct spmat_node** rows = get_private(g->groupSubmatrix);
	struct spmat_node* cur;
	printf("\n %s", "group size: ");
	printf("%d", n);
	for (i = 0; i < n; i++) {
		cur = rows[i];
		if (cur == NULL)
			printf("\n row %d is NULL", i);
		else {
			printf("\n row %d :", i);
			while (cur != NULL) {
				printf("\n data %d ", cur->data);
				printf("name %d ", cur->node_name);
				cur = cur->next;
			}
		}
	}

}

double calcModularity(struct graph* graph, double* vectorS,
		struct divisionGroup* g, double sumKiSi, double* curSAS) {
	double SBS;
	double SAS;
	double* modularity_temp;
	double* AtimesS;

	AtimesS = (double*) malloc(g->groupSize * sizeof(double));
	if (AtimesS == NULL)
		panic(ERROR_MALLOC_FAILED);

	mult_ll(g->groupSubmatrix, vectorS, AtimesS); /*A times S*/

	SAS = dotProduct(vectorS, AtimesS, g->groupSize); /*SAS*/
	modularity_temp = modularityTimesS(graph, vectorS, g, sumKiSi); /*B^ times S */

	SBS = dotProduct(vectorS, modularity_temp, g->groupSize); /* Stimes B^  S*/

	free(modularity_temp); /*free temp calc*/
	free(AtimesS); /*free temp calc*/
	*curSAS = SAS;
	return SAS - SBS;
}

void unmovedLoop(struct graph* graph, struct divisionGroup* g,
		struct node* currentNode, double* vectorS, double sumKiSi, double Q0,
		int *switchFirstUnmovedIteration, double *maxModularityChange,
		int *indexOfBiggestIncrease, struct node** prevOfBiggest) {

	/*clock_t start, end;*/

	double modChange;
	double prevSAS;
	struct node* prev;
	prev = NULL;
	/*start = clock();
	 srand(time(NULL));*/
	/*while going through all UNMVOED, they all have the same SAS. so we calculate it once here:*/
	prevSAS = calcSAS(g, vectorS);

	while (currentNode != NULL) {

		flipVectorEntry(vectorS, currentNode->data.num);

		/*modChange = calculateChangeModularity(graph, g, vectorS,
				sumKiSi, Q0, currentNode->data.num, &prevSAS);*/ /*TODO: fix or delete*/
		modChange = calculateChangeModularityWithPrevSas(graph, g, vectorS,
						sumKiSi, Q0, currentNode->data.num, prevSAS);

		if (*switchFirstUnmovedIteration == 1
				|| modChange > *maxModularityChange) {
			*maxModularityChange = modChange;
			*indexOfBiggestIncrease = currentNode->data.num; /*finding "k" of biggestIncrease;*/
			/*maxImproveScore = maxModularityChange;*/
			*prevOfBiggest = prev;
			*switchFirstUnmovedIteration = 0;
		}
		flipVectorEntry(vectorS, currentNode->data.num);

		prev = currentNode;
		currentNode = currentNode->next;
	}
	/*end = clock();*/
	/*printf("unmoved LOOP took %f seconds\n", ((double) (end - start) / CLOCKS_PER_SEC));*/

}

/*ODO: is DeltaModularity double int long?*/
void modularityMaximization(struct graph* graph, double* vectorS,
		struct divisionGroup* g) {

	double modularityChange, Q0, maxModularityChange, maxImprovedIndex = 0,
			sumKiSi, curMax, curSAS;
	int i, indexOfBiggestIncrease, switchFirstUnmovedIteration = 1;
	struct node* unmoved;
	struct node* removedFromUnmoved;
	struct node* currentNode;
	struct node* removedNode;
	struct node* prevOfBiggest;
	double* improvedVector;
	int* indiceVector;
	clock_t start, end;

	start = clock();
	srand(time(NULL));

	improvedVector = (double*) malloc(g->groupSize * sizeof(double));
	indiceVector = (int*) malloc(g->groupSize * sizeof(int));

	removedFromUnmoved = NULL;
	unmoved = createUnmovedList(g->groupSize);

	do {
		/*improving delta Q by moving ONE index*/

		sumKiSi = sumOfDegreeByVectorS(graph, vectorS, g);/* first KiSI */

		for (i = 0; i < g->groupSize; i++) {
			currentNode = unmoved;
			prevOfBiggest = NULL;
			if (i == 0)
				Q0 = calcModularity(graph, vectorS, g, sumKiSi, &curSAS);
			else
				Q0 = Q0 + calculateChangeModularity(graph, g, vectorS, sumKiSi,
								Q0, indexOfBiggestIncrease, &curSAS);
			switchFirstUnmovedIteration = 1; /*indicator that says: we need to set i = 0 as currentMax*/
			/*change in Modularity is Q1-Q0 (For Q1 with biggest modularity). we wish to find Q1 so it is Q1-Q0+Q0 ^
			 finding vertex with maximal increase in modularity*/
			sumKiSi = sumOfDegreeByVectorS(graph, vectorS, g);/* only AFTER CALC MOD we change KiSi */

			/*loop 6-10 in psuedo-code is in unmovedLoop*/
			unmovedLoop(graph, g, currentNode, vectorS, sumKiSi, Q0,
					&switchFirstUnmovedIteration, &maxModularityChange,
					&indexOfBiggestIncrease, &prevOfBiggest);
			/*moving vertex with maximal increase in modularity*/
			flipVectorEntry(vectorS, indexOfBiggestIncrease);

			unmoved = removeFromUnmoved(prevOfBiggest, unmoved,
					&removedNode);
			/*add removed node from unmoved to removeFromUnmoved list*/
			removedFromUnmoved = addToList(removedFromUnmoved,removedNode);
			/*updating vector of indice to save the index we now moved:*/
			indiceVector[i] = indexOfBiggestIncrease;
			/*updating the current 'state score'*/
			updateImprovedVector(improvedVector, i, maxModularityChange,
					&maxImprovedIndex, &curMax); /*incrementing scores*/

		}
		/*switch between lists*/
		unmoved = removedFromUnmoved;
		removedFromUnmoved = NULL;

		modularityChange = curMax;
		updateS(vectorS, indiceVector, maxImprovedIndex, g->groupSize);
		/*printf("%f \n",modularityChange);*/

		/*if all were flipped, we are the same - so loop must stop, row 26 in psuedo-code*/
		if (maxImprovedIndex == (g->groupSize) - 1)
			modularityChange = 0;

	} while (modularityChange > epsilon);
	end = clock();
	printf("ModMax took %f seconds\n",
			((double) (end - start) / CLOCKS_PER_SEC));
	free(improvedVector);
	free(indiceVector);
	free(unmoved);
}

