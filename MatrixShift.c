#include <stdio.h>
#include <stdlib.h>

/* method to calculate the 1-norm of matrix mat.
//TODO: check if double is necessary */
double one_norm(struct graph* graph, int* vectorS,
		struct divisionGroup* g)
{
	double maxColumn = 0, currentSum;
	int i;
	// iterate through columns
	for(i = 0;i<g->groupSize;i++)
	{
		currentSum = columnSum(graph,g,i);
		if (maxColumn < currentSum) maxColumn = currentSum;
	}
	return maxColumn;
}

double columnSum(struct graph* graph,
		struct divisionGroup* g,int column)
{
	double sum;
	int i;
	for(i = 0; i g->groupSize;i++)
	{
		//iterate over all rows
		sum = sum + g->groupSubmatrix.get(i,column)+g->sumOfRows[column];// TODO: yuval
		sum = sum - ((graph->vectorDegrees[i]*graph->vectorDegrees[column])/graph->M);
	}
	return sum;
}

/*given a matrix mat,and a scalar, add to mat Id*scalar */
//TODO: make sure it is okay to alter the matrix.
void AddScalarMatrix(double* mat,int rows,int cols,double scalar)
{
	int i;
	int j;
	double* matPointer;
	for (i = 0;i < rows;i++)
	{
		matPointer=mat+(i*cols)+i;
		*matPointer=*matPointer+scalar;
	}
}

void matrixShift(double* mat,int rows,int cols)
{
	double one_norm = one_norm(mat,rows,cols);
	AddScalarMatrix(mat,rows,cols,one_norm);
}
