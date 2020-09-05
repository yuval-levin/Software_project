#include<stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "Algorithm2.h"
#include "ModularityMaximization.h"

//TODO: is include file.c ok? or should we do headers?
//TODO: add checks for all mallocs.
//adds in start
void add_groupDivision(struct division* D, struct divisionGroup* g) {
	struct node* add = (struct node*) malloc(sizeof(struct node));
	if (add == NULL)
		exit(1); //TODO: print error before exit.
	add->data.group = g;
	add->next = NULL;
	if (D->len == 0)
		D->divisions = add;
	else
		add->next = D->divisions;

	D->len = (D->len) + 1;
}

/* removes  first group from D and returns it*/
struct divisionGroup* removeFirstGroup(struct division* D) {
	struct node* group = D->divisions;
	struct node* nextGroup = group->next;

	group->next = NULL;
	D->divisions = nextGroup;
	D->len = (D->len) - 1;

	return group->data.group;
}

void updateDivisionPostSplit(struct divisionGroup* g, struct division* P,
		struct division* O) {
	if (g->groupSize == 1)
		add_groupDivision(O, g);
	else
		add_groupDivision(P, g);

}
//make sure s freees g
// if there's a group size = 0, g2 = NULL, g1 = g;
void splitByS(int* vectorS, struct divisionGroup* g1, struct divisionGroup* g2) {
	//TODO
}
struct division* new_division() {
	struct division* D = (struct division*) malloc(sizeof(struct division));
	if (D == NULL)
		exit(1); //TODO: print error before exit.
	D->len = 0;
	D->divisions = NULL;
}


/*PARAMS : n is number of nodes in graph
 * DESC : Returns a divisonGroup with all nodes in graph
 */
struct divisionGroup* createTrivialDivision(int n, struct graph* inputGraph) {
	int i;
	struct divisionGroup* group = (struct divisionGroup*)malloc(sizeof(struct divisionGroup));
	if (group == NULL)
		exit(1); //TODO: print error before exit.
	group->groupSize = n;
	group->groupSubmatrix = (inputGraph->A);
	group->sumOfRows = (int*) malloc(n * sizeof(int));
	for (i = 0; i < n; i++) {
		group->sumOfRows[i] = 0;
	}
	return group;
}


struct division* Algorithm3(int numOfNodes, struct graph inputGraph) {
	struct divisionGroup* g;
	struct divisionGroup* g1;
	struct divisionGroup* g2;
	int* vectorS;
	struct division* P = new_division();
	struct division* O = new_division();
	add_groupDivision(P, createTrivialDivision(numOfNodes, &inputGraph));

	while (P->len > 0) {
		g = removeFirstGroup(P);
		vectorS = (int*) malloc(g->groupSize * sizeof(int)); //vectorS is size of group g.
		if (vectorS == NULL)
			exit(1); //TODO: print error before exit.
		Algorithm2(vectorS, g,&inputGraph);
		ModularityMaximization(&inputGraph,vectorS, g);
		splitByS(vectorS, g1, g2);

		if (g2 == NULL)
			add_groupDivision(O, g1);
		else {
			updateDivisionPostSplit(g1, P, O);
			updateDivisionPostSplit(g2, P, O);
		}
		free(vectorS);
	}
	return O;

}

