#include <stdio.h>
#include <stdlib.h>

/* method to calculate the 1-norm of matrix mat.
//TODO: check if double is necessary */
double one_norm(double *mat,int rows,int cols)
{
	int i;
	int j;
	double* matPointer;
	int maxSum = 0;
	int currentSum;
	
	for (i = 0; i < cols; i++) 
	{
		currentSum = 0 ;
		for (j = 0; j < rows; j++)
		{
			matPointer=mat+(i*cols);
			currentSum+=*matPointer;
		}
		if(currentSum > maxSum)
		{
			maxSum = currentSum;
		}	
	}
	return maxSum;
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
