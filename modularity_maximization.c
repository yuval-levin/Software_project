#include<stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "error_codes.h"
#include "spmat.h"

/* calculates the result of multiplying the changedIndex-th row of
 * mat by vectorS */
static double calc_ai_si(double* vectorS, int changedIndex, struct _spmat* mat) {
	struct spmat_node* curNode;
	curNode = get_private(mat)[changedIndex];

	return multiply_vector(curNode, vectorS);
}

double dot_product(double* vec1, double* vec2, int length) {
	/*dot product of vectors a and b*/
	int k;

	double dot = 0;

	for (k = 0; k < length; k++) {
		dot += ((*vec1) * (*vec2));
		vec1 += 1;
		vec2 += 1;
	}
	return dot;
}

/*
 * helper function:
 * receives a vector and an entry index,
 * "flips" the value at entry index from negative to positive or vice versa.
 * */
static void flip_vector_entry(double* vector, int entry) {
	vector[entry] = vector[entry] * (-1);
}

/* Helper function:
 * We wish for S be in the same state it was when we reached maxImproved Score.
 * so we reverse everything that came after it*/
static void update_s(double* vectorS, int* indiceVector, int maxImprovedIndex,
		int length) {
	int i;
	for (i = maxImprovedIndex + 1; i < length; i++) {
		flip_vector_entry(vectorS, indiceVector[i]);
	}
}

/* As seen in psuedo-code for Algo4 row 14,
 * We update "improvedVector" by Aggregating previous improvements with current improvement.
 * We also keep track of index of maximum improvement !
 * */
static void update_improved_vector(double* improvedVector, int entryIndex, double score,
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

/*helper function for handling linkedlist:
 * adds to linkedlist "list" (represented by its first node)
 * the node "node" at the beginning of list.
 * now list's first node is given "node".
 * returns start of new list.
 */
static struct node* add_to_list(struct node* list, struct node* node) {
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
static struct node* remove_from_unmoved(struct node* prevOfBiggest, struct node* unmoved,
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

/*Calculates Change in Modularity using previous SAS*/
static double calculate_change_modularity_with_prev_sas(struct graph* graph,
		struct divisionGroup* g, double* vectorS, double sumKiSi,
		double prevModularity, int changedIndex, double* previousSAS,
		int update) {
	double newModularityY;
	double currentSAS, sumAiSi;
	double* degreeDividedByM;
	int nodeNum, degree, vectorSChangedIndex;

	nodeNum = g->groupMembers[changedIndex];
	degree = graph->vectorDegrees[nodeNum];
	degreeDividedByM = graph->degreesDividedByM;
	vectorSChangedIndex = vectorS[changedIndex]; /* entry value AFTER FLIP */

	sumAiSi = calc_ai_si(vectorS, changedIndex, g->groupSubmatrix);
	currentSAS = *previousSAS + 4 * vectorSChangedIndex * sumAiSi;

	if (update > 0)/* indicator: if positive, we update the SAS value*/
	{
		*previousSAS = currentSAS; /*update SAS*/
	}
	newModularityY = prevModularity
			- 4 * vectorSChangedIndex * (degreeDividedByM[nodeNum])
					* (sumKiSi + (degree * vectorSChangedIndex)) + 0 * sumAiSi;
	newModularityY = newModularityY + (currentSAS - *previousSAS);

	return newModularityY - prevModularity;
}

/*TODO: add explanation.*/
static double* second_argument_in_calc(struct graph* graph, struct divisionGroup* g,
		double sumKiSi) {
	int i;
	double* KiDividedByMPlusSum;
	double* dividedyByM = graph->degreesDividedByM;
	KiDividedByMPlusSum = (double*) malloc(g->groupSize * sizeof(double));

	if (KiDividedByMPlusSum == NULL)
		panic(ERROR_MALLOC_FAILED);
	/*two iterations are a must, cause we need to find sum first..*/
	for (i = 0; i < g->groupSize; i++) {
		KiDividedByMPlusSum[i] = (dividedyByM[i]) * sumKiSi;
	}
	return KiDividedByMPlusSum;

}

double* modularity_times_s(struct graph* graph, double* vectorS,
		struct divisionGroup* g, double sumKiSi) {
	int i;
	double* KiDividedByMPlusSum;
	int* sumOfRows = g->sumOfRows;
	double* resVec = (double*) malloc(g->groupSize * sizeof(double));
	if (resVec == NULL)
		panic(ERROR_MALLOC_FAILED);

	KiDividedByMPlusSum = second_argument_in_calc(graph, g, sumKiSi);

	for (i = 0; i < g->groupSize; i++) {
		resVec[i] = (sumOfRows[i] * vectorS[i]) + KiDividedByMPlusSum[i];
	}

	free(KiDividedByMPlusSum);

	return resVec;
}

/*calculate kisi*/
double sum_of_degree_by_vector_s(struct graph* graph, double* vectorS,
		struct divisionGroup* g) {
	double sum = 0;
	double* vecDegrees = graph->vectorDegrees;
	int* groupMembers = g->groupMembers;
	int i;
	for (i = 0; i < g->groupSize; i++) { /*we don't know how many are really in group, since it sparse. so we use while */
		/*vectorDegrees is size Of number of nodes in the original A matrix*/
		sum = sum + (vectorS[i] * vecDegrees[groupMembers[i]]); /*vectorS is size of g, we use i*/
	}
	return sum;
}

static double sum_of_degree_by_vector_s_with_prev(struct graph* graph,
		struct divisionGroup* g, int changedIndex, double* vectorS,
		double prevKiSi) {
	double* vecDegrees = graph->vectorDegrees;
	int* groupMembers = g->groupMembers;
	double updatedKiSi = 0;

	updatedKiSi = prevKiSi
			+ 2 * vectorS[changedIndex]
					* vecDegrees[groupMembers[changedIndex]];

	return updatedKiSi;
}

static struct node* append_to_init_unmoved(struct node* prev, int index) {
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

static struct node* create_unmoved_list(int sizeOfg) {
	int i;
	struct node* head = NULL;
	struct node* prev = NULL;
	for (i = 0; i < sizeOfg; i++) {
		prev = append_to_init_unmoved(prev, i);
		if (i == 0)
			head = prev;
	}
	return head;
}

static void free_unmoved_list(struct node* unmovedListHead, int sizeOfg) {
	int i;
	struct node* nextNode;
	struct node* curNode = unmovedListHead;
	for (i = 0; i < sizeOfg; i++) {
		nextNode = curNode->next;
		free(curNode);
		curNode = nextNode;
	}
}

static double calc_modularity(struct graph* graph, double* vectorS,
		struct divisionGroup* g, double sumKiSi, double* curSAS) {
	double SBS;
	double SAS;
	double* modularity_temp;
	double* AtimesS;

	AtimesS = (double*) malloc(g->groupSize * sizeof(double));
	if (AtimesS == NULL)
		panic(ERROR_MALLOC_FAILED);

	mult_ll(g->groupSubmatrix, vectorS, AtimesS); /*A times S*/

	SAS = dot_product(vectorS, AtimesS, g->groupSize); /*SAS*/
	modularity_temp = modularity_times_s(graph, vectorS, g, sumKiSi); /*B^ times S */

	SBS = dot_product(vectorS, modularity_temp, g->groupSize); /* Stimes B^  S*/

	free(modularity_temp); /*free temp calc*/
	free(AtimesS); /*free temp calc*/
	*curSAS = SAS;
	return SAS - SBS;
}

static void unmoved_loop(struct graph* graph, struct divisionGroup* g,
		struct node* currentNode, double* vectorS, double sumKiSi, double Q0,
		int *switchFirstUnmovedIteration, double *maxModularityChange,
		int *indexOfBiggestIncrease, struct node** prevOfBiggest,
		double* curSAS) {

	/*clock_t start, end;*/

	double modChange;
	double prevSAS;
	struct node* prev;
	prev = NULL;
	/*start = clock();
	 srand(time(NULL));*/
	/*while going through all UNMVOED, they all have the same SAS. so we calculate it once here:*/
	/*prevSAS = calcSAS(g, vectorS);*/
	prevSAS = *curSAS;
	while (currentNode != NULL) {

		flip_vector_entry(vectorS, currentNode->data.num);

		modChange = calculate_change_modularity_with_prev_sas(graph, g, vectorS,
				sumKiSi, Q0, currentNode->data.num, &prevSAS, 0);

		if (*switchFirstUnmovedIteration == 1
				|| modChange > *maxModularityChange) {
			*maxModularityChange = modChange;
			*indexOfBiggestIncrease = currentNode->data.num; /*finding "k" of biggestIncrease;*/
			/*maxImproveScore = maxModularityChange;*/
			*prevOfBiggest = prev;
			*switchFirstUnmovedIteration = 0;
		}
		flip_vector_entry(vectorS, currentNode->data.num);

		prev = currentNode;
		currentNode = currentNode->next;
	}
}

/*TODO: is DeltaModularity double int long?*/
void modularity_maximization(struct graph* graph, double* vectorS,
		struct divisionGroup* g) {

	double modularityChange, Q0, maxModularityChange, maxImprovedIndex = 0,
			sumKiSi, curMax, curSAS;
	int i, indexOfBiggestIncrease,groupSize, switchFirstUnmovedIteration = 1;
	struct node* unmoved;
	struct node* removedFromUnmoved;
	struct node* currentNode;
	struct node* removedNode;
	struct node* prevOfBiggest;
	double* improvedVector;
	int* indiceVector;
	groupSize = g->groupSize;

	improvedVector = (double*) malloc(groupSize * sizeof(double));
	indiceVector = (int*) malloc(groupSize * sizeof(int));

	removedFromUnmoved = NULL;
	unmoved = create_unmoved_list(groupSize);

	do {
		/*improving delta Q by moving ONE index*/

		sumKiSi = sum_of_degree_by_vector_s(graph, vectorS, g);/* first KiSI */

		for (i = 0; i < g->groupSize; i++) {
			currentNode = unmoved;
			prevOfBiggest = NULL;
			if (i == 0)
				Q0 = calc_modularity(graph, vectorS, g, sumKiSi, &curSAS);
			else
				Q0 = Q0
						+ calculate_change_modularity_with_prev_sas(graph, g,
								vectorS, sumKiSi, Q0, indexOfBiggestIncrease,
								&curSAS, 1);
			switchFirstUnmovedIteration = 1; /*indicator that says: we need to set i = 0 as currentMax*/
			/*change in Modularity is Q1-Q0 (For Q1 with biggest modularity). we wish to find Q1 so it is Q1-Q0+Q0 ^
			 finding vertex with maximal increase in modularity*/
			/*	sumKiSi = sumOfDegreeByVectorS(graph, vectorS, g); only AFTER CALC MOD we change KiSi */
			/* if i = 0 there is no "index of biggestincrease,*/

			if (i != 0)
				sumKiSi = sum_of_degree_by_vector_s_with_prev(graph, g,
						indexOfBiggestIncrease, vectorS, sumKiSi);
			/*loop 6-10 in psuedo-code is in unmovedLoop*/
			unmoved_loop(graph, g, currentNode, vectorS, sumKiSi, Q0,
					&switchFirstUnmovedIteration, &maxModularityChange,
					&indexOfBiggestIncrease, &prevOfBiggest, &curSAS);
			/*moving vertex with maximal increase in modularity*/
			flip_vector_entry(vectorS, indexOfBiggestIncrease);

			unmoved = remove_from_unmoved(prevOfBiggest, unmoved, &removedNode);
			/*add removed node from unmoved to remove_from_unmoved list*/
			removedFromUnmoved = add_to_list(removedFromUnmoved, removedNode);
			/*updating vector of indice to save the index we now moved:*/
			indiceVector[i] = indexOfBiggestIncrease;
			/*updating the current 'state score'*/
			update_improved_vector(improvedVector, i, maxModularityChange,
					&maxImprovedIndex, &curMax); /*incrementing scores*/

		}
		/*switch between lists*/
		unmoved = removedFromUnmoved;
		removedFromUnmoved = NULL;

		modularityChange = curMax;
		update_s(vectorS, indiceVector, maxImprovedIndex, groupSize);/*psuedo code row 22-24*/

		/*if all were flipped, we are the same - so loop must stop, row 26 in psuedo-code*/
		if (maxImprovedIndex == (groupSize) - 1)
			modularityChange = 0;

	} while (modularityChange > epsilon);

	free(improvedVector);
	free(indiceVector);
	free_unmoved_list(unmoved,groupSize);
}

