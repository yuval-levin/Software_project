/*
 * spmat.c
 * how to create a sparse matrix:
 * spmat* mat;
 * mat = spmat_allocate_list(n);
 * create_matrix(mat, input);
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spmat.h"
#include <time.h>
#include <math.h>
#include <string.h>

/*
 * spmat_node structure which implements the linked list.
 * index contains the column index.
 */
typedef struct spmat_node {
	double data;
	int index;
	struct spmat_node *next;
} spmat_node;

void add_row_ll(struct _spmat *A, const double *row, int i);
void free_ll(struct _spmat *A);
void mult_ll(const struct _spmat *A, const double *v, double *result);


/*
 * allocating memory for rows.
 * populate the functions with linked list implementation.
 */
spmat* spmat_allocate_list(int n){
	spmat* mat;
	spmat_node** rows;

	/*allocating memory*/
	mat = (spmat*)malloc(sizeof(spmat));
	assert(mat!= NULL);
	rows = (spmat_node**)malloc(n * sizeof(spmat_node*));
	assert(rows!= NULL);

	/*initializing*/
	mat->n = n;
	mat->add_row = &add_row_ll;
	mat->free = &free_ll;
	mat->mult = &mult_ll;
	mat->private = rows;

	return mat;
}

/*
 * add_row implementation for linked list
 */
void add_row_ll(struct _spmat *A, const double *row, int i){
	add_row_of_size_n(A->n);
}

/*
 * helper for add_row_ll
 */
void add_row_of_size_n(struct _spmat *A, const double *row, int i, int n){
	int first, j;
	spmat_node *cur, *prev;
	first = 0;

	for (j = 0; j < n; j++){
		if (row[j] != 0){
			cur = (spmat_node*)malloc(sizeof(spmat_node));
			assert(cur != NULL);
			cur->data = row[j];
			cur->index = j;		/*column index*/
			cur->next = NULL;

			/*first spmat_node*/
			if (first == 0){
				((spmat_node**)A->private)[i] = cur;	/*inserted as the 1st node of the ith row*/
				first = 1;							/*no updating prev.next*/
			}
			else{
				prev->next = (spmat_node*)cur;		/*cur is not the 1st, update prev.next*/
			}
			prev = cur;
		}
	}
}

/*
 * helper for free_ll, frees the nodes of the given row recursively
 */
void free_row_ll(spmat_node* row){
	if (row != NULL){
		free_row_ll((spmat_node*)row->next);
		free(row);
	}
}

/*
 * free implementation for linked list
 */
void free_ll(struct _spmat *A){
	int i;
	spmat_node** rows;
	rows = A->private;

	for(i = 0; i< A->n; i++){
		free_row_ll(rows[i]);
	}
	free(A->private);
	free(A);
}

/*
 * mult implementation for linked list
 */
void mult_ll(const struct _spmat *A, const double *v, double *result){
	int row_ind;
	int index;
	double sum;
	spmat_node* cur_node;
	spmat_node** rows = (spmat_node** )A->private;

	/*calculate the result of each row*/
	for(row_ind = 0; row_ind < A->n; row_ind++){
		sum = 0;
		cur_node = rows[row_ind];
		
		/*sum all the relevant multiplications*/
		while(cur_node != NULL){
			index = cur_node->index;
			sum += (cur_node->data) * (v[index]);
			cur_node = (spmat_node*)cur_node->next;
		}
		result[row_ind] = sum;
	}
}

/*
 * creating the matrix one row at a time
 */
void create_matrix(struct _spmat *A, FILE* input){
	double *row;
	int i, k, n;

	n = A->n;
	fseek(input, 8, SEEK_SET); 	/*skipping dimensions*/
	row = (double*)malloc(n*sizeof(double));
	assert(row != NULL);
	/*add each row*/
	for (i = 0; i < n; i++){
		k = fread(row, sizeof(double),n, input);
		assert (k == n);
		A->add_row(A, row, i);
	}

	free(row);
	rewind(input);
}