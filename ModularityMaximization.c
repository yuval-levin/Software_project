#include<stdio.h>
#include <stdlib.h>
#include "modules.h"

//TODO: is DeltaModularity double int long?
void modularityMaximization(int* vectorS, struct divisionGroup* g)
{
	const double epsilon = 0.00001;
	double modularityChange;
	do
	{

	}while(modularityChange > epsilon);
}

struct node* createUnmovedList(int sizeOfg)
{
	int i;
	struct node* head,prev;
	head =(struct node*)malloc(sizeof(struct node*))
	head->data.num = 0;
	head->next = NULL;
	prev = head;
	for(i = 1;i < sizeOfg;i++)
}
