#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "error_codes.h"
#include "spmat.h"

double multiply_vector(struct spmat_node* curNode, const double *v) {
    double sum = 0;
    for (; curNode != NULL; curNode = curNode->next) {
        sum += v[curNode->index];
    }
    return sum;
}

/*
 * mult implementation for linked list
 */
void mult_ll(const struct _spmat *A, const double *v, double *result) {
	int rowInd;
	double sum;
	struct spmat_node* cur_node;
	struct spmat_node** rows = (struct spmat_node**) A->private;

	/*calculate the result of each row*/
	for (rowInd = 0; rowInd < A->n; rowInd++) {
		cur_node = rows[rowInd];

		/*sum all the relevant multiplications*/
		sum = multiply_vector(cur_node, v);
		result[rowInd] = sum;
	}
}

void add_row_of_size_n(struct _spmat *A, const double *row, int i, int n) {
	int j;
	struct spmat_node *cur, *prev;
	prev = NULL;

	if (n == 0) {
		((struct spmat_node**) A->private)[i] = NULL;
	}

	for (j = 0; j < n; j++) {
		cur = (struct spmat_node*) malloc(sizeof(struct spmat_node));
		if (cur == NULL)
			panic(ERROR_MALLOC_FAILED);

		cur->data = 1;
		cur->index = row[j];
		cur->node_name = row[j];
		cur->next = NULL;

		/*first struct spmat_node*/
		if (prev == NULL) {
			((struct spmat_node**) A->private)[i] = cur; /*inserted as the 1st node of the ith row*/
			/*no updating prev.next*/
		} else {
			prev->next = (struct spmat_node*) cur; /*cur is not the 1st, update prev.next*/
		}
		prev = cur;
	}
}

/*
 * add_row implementation for linked list
 */
void add_row_ll(struct _spmat *A, const double *row, int i) {
	add_row_of_size_n(A, row, i, A->n);
}

struct spmat_node** get_private(struct _spmat* mat) {
	return mat->private;
}

void set_private(struct _spmat* mat, struct spmat_node** rows) {
	mat->private = rows;
}

/*
 * helper for free_ll, frees the nodes of the given row recursively
 */
void free_row_ll(struct spmat_node* row) {
	if (row != NULL) {
		free_row_ll((struct spmat_node*) row->next);
		free(row);
	}
}

/*
 * free implementation for linked list
 */
void free_A(struct _spmat *A) {
	free(A->private);
	free(A);
}

/*
 * free implementation for linked list
 */
void free_ll(struct _spmat *A) {
	int i;
	struct spmat_node** rows;
	rows = A->private;

	if (rows != NULL) {
		for (i = 0; i < A->n; i++) {
			free_row_ll(rows[i]);
		}
	}
	free_A(A);
}

/*
 * allocating memory for rows.
 * populate the functions with linked list implementation.
 */
spmat* spmat_allocate_list(int n) {
	spmat* mat;
	struct spmat_node** rows;

	/*allocating memory*/
	mat = (spmat*) malloc(sizeof(spmat));
	if (mat == NULL)
		panic(ERROR_MALLOC_FAILED);

	rows = (struct spmat_node**) malloc(n * sizeof(struct spmat_node*));
	if (rows == NULL)
		panic(ERROR_MALLOC_FAILED);

	/*initializing*/
	mat->n = n;
	mat->add_row = &add_row_ll;
	mat->free = &free_ll;
	mat->mult = &mult_ll;
	mat->private = rows;

	return mat;
}

spmat* spmat_allocate_list_without_rows(int n) {
	spmat* mat;

	/*allocating memory*/
	mat = (spmat*) malloc(sizeof(spmat));
	if (mat == NULL)
		panic(ERROR_MALLOC_FAILED);

	/*initializing*/
	mat->n = n;
	mat->add_row = &add_row_ll;
	mat->free = &free_ll;
	mat->mult = &mult_ll;

	return mat;
}

