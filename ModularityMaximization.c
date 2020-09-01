#include<stdio.h>
#include <stdlib.h>
#include "modules.h"

//TODO: is DeltaModularity double int long?
void modularityMaximization(struct graph* graph,int* vectorS, struct divisionGroup* g)
{
	const double epsilon = 0.00001;
	double modularityChange,Q0,Q1,maxModularityChange,maxImprovedIndex = 0,maxImproveScore;
	int i, indexOfBiggestIncrease;
	struct node* unmoved,currentNode,prev,prevOfBiggest;
	double* improvedVector = (double*)malloc(g->groupSize*sizeof(double));
	int* indiceVector = (int*)malloc(g->groupSize*sizeof(int));

	do
	{
		//improving delta Q by moving ONE index
		unmoved = createUnmovedList(g->groupSize);

		for(i = 0;i< g->groupSize;i++)
		{
			currentNode = unmoved;
			prev = NULL, prevOfBiggest = NULL;
			Q0 = dotProduct(vectorS,modularityTimesS(graph,vectorS,g));
			//finding vertex with maximal increase in modularity
			while (currentNode != NULL)
			{
				flipVectorEntry(vectorS,currentNode->data.num);
				Q1 = calculateChangeModularity(); //helper function in O1, todo

				if (i == 0 || Q1 > maxModularityChange)
					{
						maxModularityChange = Q1 , indexOfBiggestIncrease = i;
						maxImproveScore = maxModularityChange;
						prevOfBiggest = prev;
					}
				flipVectorEntry(vectorS,currentNode->data.num);

				prev = currentNode;
				currentNode = currentNode->next;
			}
			//moving vertex with maximal increase in modularity
			flipVectorEntry(vectorS,indexOfBiggestIncrease);
			removeFromUnmoved(prevOfBiggest,unmoved);
			//updating vector of indice to save the index we now moved:
			indiceVector[i]=indexOfBiggestIncrease;
			//updating the current 'state score'
			updateImprovedVector(improvedVector,i,maxModularityChange); // incrementing scores
			if (improvedVector[i] > maxImproveScore) maxImproveScore = improvedVector[i], maxImprovedIndex=i;
		}
		modularityChange = maxImproveScore;
		updateS(vectorS,indiceVector,maxImprovedIndex,g->groupSize);

	}while(modularityChange > epsilon);

	free(improvedVector);
	free(indiceVector);
	free(unmoved);
}

/* we wish for S be in the same state it was when we reached maxImproved Score.
 * so we reverse everything that came after it*/
void updateS(int* vectorS,int* indiceVector,int maxImprovedIndex,int length)
{
	int i;
	for(i = maxImprovedIndex +1 ;i < length; i++)
	{
		flipVectorEntry(vectorS,indiceVector[i]);
	}
}
void flipVectorEntry(int* vector, int entry)
{
	vector[entry] = vector[entry]*(-1);
}

void updateImprovedVector(int* improvedVector, int entryIndex,double score)
{
	if (entryIndex == 0) improvedVector[0] = score;
	else improvedVector[entryIndex] = improvedVector[entryIndex-1]+score;
}

void removeFromUnmoved(struct node* prevOfBiggest, struct node* unmoved)
{
	struct node* removedNode;
	if (prevOfBiggest == NULL) unmoved = unmoved->next , removedNode = unmoved; //todo: make sure this works
	else prevOfBiggest->next = removedNode = prevOfBiggest->next, prevOfBiggest->next->next;
	free(removedNode);
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
