#include<stdio.h>
#include <stdlib.h>
#include "modules.h"

//TODO: is DeltaModularity double int long?
void modularityMaximization(int* vectorS, struct divisionGroup* g)
{
	const double epsilon = 0.00001;
	double modularityChange,Q0;
	int i;
	struct node* unmoved = createUnmovedList(g->groupSize);

	do
	{
		for(i = 0;i< g->groupSize;i++)
		{

		}

	}while(modularityChange > epsilon);
}

struct node* createUnmovedList(int sizeOfg)
{
	int i;
	struct node* head,prev = NULL;
	for(i = 0;i < sizeOfg;i++)
		{
			prev = appendToList(prev,i);
			if (i == 0) head = prev;
		}
	return head;
}

struct node* appendToList(struct node* prev, int index)
{
	struct node* current;

	current = (struct node*)malloc(sizeof(struct node*));
	current->data.num = index;
	current->next = NULL;
	if (prev != NULL) prev->next = current;

	return current;
}
