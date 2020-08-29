/*
 * spmat.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spmat.h"
#include <time.h>
#include <math.h>
#include <string.h>

/*
 * node structure which implements the linked list.
 * index contains the column index.
 */
typedef struct spmat_node {
	double data;
	int index;
	struct node *next;
} spmat_node;

/*
 * linked list implementation of spmat
 */
typedef struct spmat_ll {
	node** rows; 	/*a list of nodes, each node is the first node
					in the corresponding row*/
} spmat_ll;

/*
 * arrays implementation of spmat
 */
typedef struct spmat_ar {
	double* values;		/*compressed array of the non-zero values*/
	int* colind;		/*colind[i] = column of values[i]
	 	 	 	 	 	 in the original matrix*/
	int* rowptr;		/*rowptr[i] = index in the compressed array of
						the first non-zero value from the ith row onwards*/
	int val_index;		/*index of the first non-updated cell in values*/
	int nnz;			/*number of non-zero values*/
} spmat_ar;

void add_row_ll(struct _spmat *A, const double *row, int i);
void free_ll(struct _spmat *A);
void mult_ll(const struct _spmat *A, const double *v, double *result);
void add_row_ar(struct _spmat *A, const double *row, int i);
void free_ar(struct _spmat *A);
void mult_ar(const struct _spmat *A, const double *v, double *result);
void rowptr_update(struct spmat_ar *A, int n);
double mult_row_by_val(struct spmat_ar *mat_ar, int i_start, int i_end, const double *v);

/*
 * allocating memory for rows.
 * populate the functions with linked list implementation.
 */
spmat* spmat_allocate_list(int n){
	spmat *mat;
	spmat_ll *mat_ll;

	/*allocating memory*/
	mat = (spmat*)malloc(sizeof(spmat));
	assert(mat!= NULL);
	mat_ll = (spmat_ll*)malloc(sizeof(spmat_ll));
	assert(mat_ll!= NULL);

	mat_ll->rows = (node**)malloc(n * sizeof(node*));
	assert(mat_ll->rows != NULL);

	/*initializing*/
	mat->n = n;
	mat->add_row = add_row_ll;
	mat->free = free_ll;
	mat->mult = mult_ll;
	mat->private = &mat_ll;

	return mat;
}

/*
 * add_row implementation for linked list
 */
void add_row_ll(struct _spmat *A, const double *row, int i){
	int n;
	int first, j;
	node *cur, *prev;
	spmat_ll* mat_ll;

	n = A->n;
	first = 0;
	mat_ll = A->private;

	for (j = 0; j < n; j++){
		if (row[j] != 0){
			cur = (node*)malloc(sizeof(node));
			assert(cur != NULL);
			cur->data = row[j];
			cur->index = j;		/*column index*/
			cur->next = NULL;

			/*first node*/
			if (first == 0){
				((node**)mat_ll->rows)[i] = cur;	/*inserted as the 1st node of the ith row*/
				first = 1;							/*no updating prev.next*/
			}
			else{
				prev->next = (node*)cur;		/*cur is not the 1st, update prev.next*/
			}
			prev = cur;
		}
	}
}

/*
 * free implementation for linked list
 */
void free_ll(struct _spmat *A){
	int n;
	node *cur, *next, *p_rows;
	node **rows;
	spmat_ll *mat_ll;

	mat_ll = A->private;
	rows = mat_ll->rows;
	n = A->n;

	/*free all nodes in rows*/
	for (p_rows = rows[0]; p_rows <= rows[n-1]; p_rows++){
		cur = p_rows;
		/*free nodes in each row*/
		while (cur != NULL){
			next = cur->next;
			free(cur);
			cur = next;
		}
	}

	/*free rows array*/
	free(rows);
	free(mat_ll);
	free(A);
}

/*
 * mult implementation for linked list
 */
void mult_ll(const struct _spmat *A, const double *v, double *result){
	int n, iter;
	int index;
	double sum;
	double *p_res;
	node *cur_node, *cur_row;
	node **rows;
	spmat_ll *mat_ll;

	mat_ll = A->private;
	rows = (node**)mat_ll->rows;
	n = A->n; 
	p_res = result;
	printf("TEST mult 1\n");
	iter = 0;
	printf(cur_row);
	cur_row = *rows;
	printf("TEST mult 2\n");

	/*calculate the result of each row*/
	for (cur_row = rows[0]; cur_row <= rows[n-1]; cur_row++){
		sum = 0;
		cur_node = cur_row;
		printf("TEST mult 3\n");

		/*sum all the relevant multiplications*/
		while (cur_node != NULL){
			index = cur_node->index;
			sum += cur_node->data * (v[index]);
			cur_node = (node*)cur_node->next;
			printf("TEST mult iter %d\n", iter);
		}

		/*sum is updated, place it in result and increment p_res*/
		*p_res = sum;
		p_res++;
	}
}

/*
 * allocating memory for the arrays.
 * populate the functions with arrays implementation.
 */
spmat* spmat_allocate_array(int n, int nnz){
	spmat *mat;
	spmat_ar *mat_ar;

	/*allocating memory*/
	mat = (spmat*)malloc(sizeof(spmat));
	assert(mat!= NULL);
	mat_ar = (spmat_ar*)malloc(sizeof(spmat_ar));
	assert(mat_ar!= NULL);

	mat_ar->values = (double*)malloc(nnz * sizeof(double));
	assert(mat_ar->values != NULL);
	mat_ar->colind = (int*)malloc(nnz * sizeof(int));
	assert(mat_ar->colind != NULL);
	mat_ar->rowptr = (int*)malloc((n+1) * sizeof(int));
	assert(mat_ar->rowptr != NULL);
	mat_ar->val_index = 0;
	mat_ar->nnz = nnz;

	/*initializing*/
	mat->n = n;
	mat->add_row = add_row_ar;
	mat->free = free_ar;
	mat->mult = mult_ar;
	mat->private = mat_ar;

	return mat;
}

/*
 * add_row implementation for arrays
 */
void add_row_ar(struct _spmat *A, const double *row, int i){
	int n, k, cur_row_ptr, j, h;
	double *values;
	int *colind, *rowptr;
	spmat_ar* mat_ar;

	n = A->n;
	mat_ar = A->private;
	k = mat_ar->val_index;		/*index of the first non-updated cell in values*/
	values = mat_ar->values;
	colind = mat_ar->colind;
	rowptr = mat_ar->rowptr;

	values += k;				/*points to the first cell to be updated*/
	colind +=k;
	cur_row_ptr = -1;			/*index of the first non-zero in the row*/
	j = 0;

	for (h = 0; h < n; h++){
		if (row[h] != 0){
			*values = row[h];
			*colind = j;		/*column num*/

			/* update cur_row_ptr once,
			 * the first non-zero val in the row*/
			if (cur_row_ptr == -1){
				cur_row_ptr = k;
			}
			values++;
			colind++;
			k++;
		}
		j++;
	}

	rowptr[i] = cur_row_ptr;	/*if the row contains only 0, put -1*/
	mat_ar->val_index = k;		/*next cell to be updated is k*/

	/*last row, handle rowptr*/
	if (i == n-1 && n != 0){
		rowptr_update(mat_ar, n);
	}
}

/*
 * helper for add_row_ar,
 * updates rowptr.
 */
void rowptr_update(struct spmat_ar *A, int n){
	int *p, *rowptr;

	rowptr = A->rowptr;
	rowptr[n] = A->nnz;

	/* if rowptr[i] == -1, the ith row contains only zeros.
	 * its value in rowptr should be rowptr[i+1]*/
	for (p = rowptr + n-1; p >= rowptr; p--){
		if (*p == -1){
			*p = *(p+1);
		}
	}
}

/*
 * free implementation for arrays
 */
void free_ar(struct _spmat *A){
	double *values;
	int *colind, *rowptr;
	spmat_ar *mat_ar;

	mat_ar = A->private;
	values = mat_ar->values;
	colind = mat_ar->colind;
	rowptr = mat_ar->rowptr;

	free(values);
	free(colind);
	free(rowptr);
	free(mat_ar);
	free(A);
}

/*
 * mult implementation for arrays
 */
void mult_ar(const struct _spmat *A, const double *v, double *result){
	int n, i;
	int *rowptr, *p_rowptr;
	spmat_ar* mat_ar;

	n = A->n;
	mat_ar = A->private;
	rowptr = mat_ar->rowptr;

	p_rowptr = rowptr + n-1;		/*points to the last row*/
	result = result + n-1;			/*we will calc result from last entry to first*/

	/*for each row (from last to first), calc res*/
	for (i = n-1; i >= 0; i--){

		/*if the current row contains only 0, res = 0*/
		if(*p_rowptr == *(p_rowptr+1)){
			*result = 0;
		}
		else{
			*result = mult_row_by_val(mat_ar, *p_rowptr, *(p_rowptr+1), v);
		}

		result--;
		p_rowptr--;
	}
}

/*
 * helper for mult_ar,
 * calculates a specific entry in result.
 */
double mult_row_by_val(struct spmat_ar *mat_ar, int i_start, int i_end, const double *v){
	int k, j;
	double sum;
	double *values;
	int *colind;

	values = mat_ar->values;
	colind = mat_ar->colind;
	values = &values[i_start];
	colind = &colind[i_start];
	sum = 0;

	for (k = i_start; k < i_end; k++){
		j = *colind;					/*the relevant entry in v is the index colind[k]*/
		sum += (*values) * (v[j]);		/*mult values[k] by the relevant entry in v*/
		values++;
		colind++;
	}
	return sum;
}
