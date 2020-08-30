#include<stdio.h>
#include <stdlib.h>
#include "modules.h"

//TODO: is DeltaModularity double int long?
void modularityMaximization(struct graph* graph,int* vectorS, struct divisionGroup* g)
{
	const double epsilon = 0.00001;
	double modularityChange,Q0,Q1,maxModularityChange;
	int i, indexOfBiggestIncrease;
	struct node* unmoved,currentNode,prevOfBiggest = NULL;

	do
	{
		unmoved = createUnmovedList(g->groupSize);
		currentNode = unmoved;
		for(i = 0;i< g->groupSize;i++)
		{
			Q0 = dotProduct(vectorS,modularityTimesS(graph,vectorS,g));
			while (currentNode != NULL)
			{
				vectorS[currentNode->data.num] = (-1)*vectorS[currentNode->data.num];
				Q1 = calculateChangeModularity(); //helper function in O1, todo
				if (i == 0 || Q1 > maxModularityChange) maxModularityChange = Q1 , indexOfBiggestIncrease = i;
				vectorS[currentNode->data.num] = (-1)*vectorS[currentNode->data.num];
				currentNode = currentNode->next;
			}

		}

	}while(modularityChange > epsilon);

	free(unmoved);
}
double calculateChangeModularity()
{
	//TODO: implement me;
	return 0;
}

double* modularityTimesS(struct graph* graph,int* vectorS, struct divisionGroup* g)
{
	int i;
	double* resVec = (double*)malloc(g->groupSize*sizeof(double));
	if(resVec == NULL) exit(1); //TODO: print error before exit.
	double* KiDividedByMPlusSum = secondArgumentInCalc(graph,vectorS,g);

	for(i = 0 ; i < g->groupSize;i++)
	{
		resVec[i] = g->sumOfRows[i]*vectorS[i] + (graph->vectorDegrees) + KiDividedByMPlusSum[i];
	}

	free(KiDividedByMPlusSum);

	return resVec;
}

//TODO: add explanation.
double* secondArgumentInCalc(struct graph* graph,int* vectorS,struct divisionGroup* g)
{
	int i,index;
	double  sum = 0;
	double* KiDividedByMPlusSum = (double*)malloc(g->groupSize * sizeof(double));
	if(KiDividedByMPlusSum == NULL) exit(1); //TODO: print error before exit.
	struct spmat_node* current = g->groupSubmatrix->private[0];

	for(i = 0;i < g->groupSize;i++)
	{
		//vectorDegrees is size Of number of nodes in the original A matrix
		sum = sum + (vectorS[i]*graph->vectorDegrees[current->index]); //vectorS is size of g, we use i
		current = current->next;
		KiDividedByMPlusSum[i] = graph->vectorDegrees[current->index] / graph->M;
	}
	//two iterations are a must, cause we need to find sum first..
	for(i = 0;i < g->groupSize;i++)
	{
		KiDividedByMPlusSum[i] = (graph->vectorDegrees[current->index] / graph->M)*sum;
	}
	return KiDividedByMPlusSum;


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
	if(current == NULL) exit(1); //TODO: print error before exit.
	current->data.num = index;
	current->next = NULL;
	if (prev != NULL) prev->next = current;

	return current;
}


double dotProduct(double* a,double* b,int col)
{
	/*dot product of vectors a and b*/
	int k;
	double* vec1 = a,vec2 = b;
	double dot = 0;

	for(k=0;k<col;k++)
	{
		dot+=((*vec1)*(*vec2));
		vec1+=1;
		vec2+=1;
	}
	return dot;
}
