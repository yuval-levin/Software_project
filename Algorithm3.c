#include<stdio.h>
#include <stdlib.h>
#include "modules.h"
//TODO: add checks for all mallocs.

struct division* Algorithm3(int numOfNodes,struct graph inputGraph)
{
	struct divisionGroup* g = NULL ,g1 = NULL,g2 = NULL;
	int* vectorS = (int*)malloc(numOfNodes*sizeof(int));

	struct division* P = new_division();
	struct division* O = new_division();
	add_groupDivision(P,createTrivialDivision(numOfNodes,&inputGraph));

	while(P->len > 0)
	{
		g = removeFirstGroup(P);
		Algorithm2(vectorS,g);
		ModularityMaximization(vectorS,g);
		splitByS(vectorS,g1,g2);

		if (g2 == NULL) add_groupDivision(O,g1);
		else
			{
				updateDivisionPostSplit(g1,P,O);
				updateDivisionPostSplit(g2,P,O);
			}
	}
	return O;

}

void updateDivisionPostSplit(struct divisionGroup* g,struct division* P,struct division* O)
{
	if (g->groupSize == 1) add_groupDivision(O,g);
	else  add_groupDivision(P,g);

}
//make sure s freees g
// if there's a group size = 0, g2 = NULL, g1 = g;
void splitByS(int* vectorS, struct divisionGroup* g1, struct divisionGroup* g2)
{
	//TODO
}
struct division* new_division()
{
	struct division D =(struct division)malloc(sizeof(struct division));
	D->len = 0;
	D->divisions = NULL;
}

//adds in start
void add_groupDivision(struct division* D,struct divisionGroup* g)
{
	struct node* add = (struct node)malloc(sizeof(struct node));
	add->data = g;
	add->next = NULL;
	if(D->len == 0) D->divisions = add;
	else add->next = D->divisions;

	D->len = (D->len) +1;
}

struct divisionGroup* removeFirstGroup(struct division* D)
{
	struct divisonGroup* group = D->divisions;
	struct divisonGroup* nextGroup = group->next;
	group->next = NULL;
	D->divisions = nextGroup;
	D->len = (D->len)-1;

	return group;
}
/*PARAMS : n is number of nodes in graph
 * DESC : Returns a divisonGroup with all nodes in graph
 */
struct divisionGroup* createTrivialDivision(int n, struct graph* inputGraph)
{
	int i;
	struct divisonGroup* group = (struct divisionGroup)malloc(sizeof(struct divisionGroup));
	group->size = n;
	group->groupSubmatrix = &(inputGraph->A);
	group->sumOfRows = (int*)malloc(n*sizeof(int));
	for (i = 0; i < n; i++)
	{
		group->sumOfRows[i] = 0;
	}
	return group;
}



