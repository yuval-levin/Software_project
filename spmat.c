/*
 * spmat.c
 */
#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include <time.h>
#include <math.h>
#include <string.h>

void add_row_ll(struct _spmat *A, const int *row, int i);
void free_ll(struct _spmat *A);
void mult_ll(const struct _spmat *A, const int *v, int *result);


/*
 * allocating memory for rows.
 * populate the functions with linked list implementation.
 */
spmat* spmat_allocate_list(int n){
	spmat* mat;
	struct spmat_node** rows;

	/*allocating memory*/
	mat = (spmat*)malloc(sizeof(spmat));
	assert(mat!= NULL);
	rows = (struct spmat_node**)malloc(n * sizeof(struct spmat_node*));
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
void add_row_ll(struct _spmat *A, const int *row, int i){
	add_row_of_size_n(A, row, i, A->n, 0);
}

struct spmat_node** get_private(struct _spmat* mat)
{
	return mat->private;
}

void set_private(struct _spmat* mat, struct spmat_node** rows){
	mat->private = rows;
}

/*
 * helper for add_row_ll
 */
void add_row_of_size_n(struct _spmat *A, const int *row, int i, int n, int is_adjacency_mat){
	int first, j;
	struct spmat_node *cur;
	struct spmat_node *prev;
	first = 0;

	for (j = 0; j < n; j++){
		if (row[j] != 0){
			cur = (struct spmat_node*)malloc(sizeof(struct spmat_node));
			assert(cur != NULL);
			if (is_adjacency_mat == 1)
			{
				cur->data = 1;
				cur->index = row[j];
			} 
			else {
				cur->data = row[j];
				cur->index = j;		/*column index*/
			}
			cur->next = NULL;

			/*first struct spmat_node*/
			if (first == 0){
				((struct spmat_node**)A->private)[i] = cur;	/*inserted as the 1st node of the ith row*/
				first = 1;							/*no updating prev.next*/
			}
			else{
				prev->next = (struct spmat_node*)cur;		/*cur is not the 1st, update prev.next*/
			}
			prev = cur;
		}
	}
}

/*
 * helper for free_ll, frees the nodes of the given row recursively
 */
void free_row_ll(struct spmat_node* row){
	if (row != NULL){
		free_row_ll((struct spmat_node*)row->next);
		free(row);
	}
}

/*
 * free implementation for linked list
 */
void free_ll(struct _spmat *A){
	int i;
	struct spmat_node** rows;
	rows = A->private;

	for(i = 0; i< A->n; i++){
		free_row_ll(rows[i]);
	}
	free(A->private);
	free(A);
}


int sumTimesVectorHelper(struct spmat_node* cur_node, int *v)
{
	int index;
	int sum = 0;
	while(cur_node != NULL){
				index = cur_node->index;
				sum += (cur_node->data) * (v[index]);
				cur_node = (struct spmat_node*)cur_node->next;
			}
	return sum;
}

/*
 * mult implementation for linked list
 */
void mult_ll(const struct _spmat *A, const int *v, int *result){
	int row_ind;
	int sum;
	struct spmat_node* cur_node;
	struct spmat_node** rows = (struct spmat_node** )A->private;

	/*calculate the result of each row*/
	for(row_ind = 0; row_ind < A->n; row_ind++){
		cur_node = rows[row_ind];
		
		/*sum all the relevant multiplications*/
		sum = sumTimesVectorHelper(cur_node, v);
		result[row_ind] = sum;
	}
}


/*
 * creating the matrix one row at a time
 */
void create_matrix(struct _spmat *A, FILE* input){
	int *row;
	int i, k, n;

	n = A->n;
	fseek(input, 8, SEEK_SET); 	/*skipping dimensions*/
	row = (int*)malloc(n*sizeof(int));
	assert(row != NULL);
	/*add each row*/
	for (i = 0; i < n; i++){
		k = fread(row, sizeof(int),n, input);
		assert (k == n);
		A->add_row(A, row, i);
	}

	free(row);
	rewind(input);
}
